using JuMP
using HiGHS
using CSV
using DataFrames
using Interpolations

# ==========================================
# 1. DATA LOADING & ALIGNMENT
# ==========================================

# --- Load 30-min datacentre data ---
df_30 = CSV.read("C:/Users/naeem/Downloads/datacentre_nsw_merged5.csv", DataFrame)

println("Datacentre CSV columns: ", names(df_30))

# Convert demand MW -> kW
demand = df_30.DATACENTRE_DEMAND .* 1000

T = length(demand)
dt = 0.5   # 30-minute timestep (hours)

# --- Load hourly solar data ---
df_solar_hourly = CSV.read(
    "C:/Users/naeem/Downloads/ninja_pv_-31.3665_144.8481_uncorrected.csv",
    DataFrame;
    header=4
)

println("Solar CSV columns: ", names(df_solar_hourly))

solar_hourly = df_solar_hourly.electricity

# Interpolate hourly -> 30-min
h_grid = 1:length(solar_hourly)
itp = linear_interpolation(h_grid, solar_hourly, extrapolation_bc=Flat())

solar_30min = [itp(t) for t in range(1, stop=length(solar_hourly), length=T)]

# ==========================================
# 2. ECONOMIC PARAMETERS (Annualized)
# ==========================================

discount_rate = 0.06

fom_pv = 21.001
fom_batt = 38.78
vom_batt = 0.0
vom_pv = 0.0

crf(i, n) = (i * (1 + i)^n) / ((1 + i)^n - 1)

k_pv_ann   = 1491.64 * crf(discount_rate, 30)   # $/kW/year
k_batt_ann = 1711.40 * crf(discount_rate, 30)   # assumed $/kWh/year

# ==========================================
# 3. OPTIMIZATION MODEL (OFF-GRID)
#    LIGHTER ALTERNATIVE:
#    - one signed battery power variable
#    - no simultaneous charge/discharge possible
#    - simplified battery efficiency treatment
# ==========================================

model = Model(HiGHS.Optimizer)
set_string_names_on_creation(model, false)

# --- Capacity variables ---
@variable(model, C_pv >= 0)      # kW
@variable(model, C_batt >= 0)    # kWh

# --- Dispatch variables ---
# p_batt[t] > 0  => battery discharging
# p_batt[t] < 0  => battery charging
@variable(model, p_batt[1:T])
@variable(model, e_batt[1:T] >= 0)
@variable(model, l_actual[1:T] >= 0)

# OBJECTIVE
@objective(model, Min,
    # PV costs
    (C_pv * k_pv_ann) +
    (C_pv * fom_pv) +

    # Battery costs (energy-based)
    (C_batt * k_batt_ann) +
    (C_batt * fom_batt)
)

alpha = 0.10

for t in 1:T

    # Power balance (OFF-GRID)
    # PV + battery net discharge = datacentre load
    @constraint(model,
        (C_pv * solar_30min[t]) + p_batt[t] == l_actual[t]
    )

    # Battery dynamics
    # Simplified signed-power model:
    #   p_batt > 0  -> discharge -> energy decreases
    #   p_batt < 0  -> charge    -> energy increases
    if t == 1
        @constraint(model,
            e_batt[t] == e_batt[T] - p_batt[t] * dt
        )
    else
        @constraint(model,
            e_batt[t] == e_batt[t-1] - p_batt[t] * dt
        )
    end

    # Battery energy limit
    @constraint(model, e_batt[t] <= C_batt)

    # 4-hour battery power limit
    @constraint(model, p_batt[t] <=  0.25 * C_batt)
    @constraint(model, p_batt[t] >= -0.25 * C_batt)

    # Load flexibility
    @constraint(model, l_actual[t] >= (1 - alpha) * demand[t])
    @constraint(model, l_actual[t] <= (1 + alpha) * demand[t])
end

# Daily conservation of datacentre energy
for d in 0:(Int(T / 48) - 1)
    day_indices = (d * 48 + 1):(d * 48 + 48)
    @constraint(model,
        sum(l_actual[t] for t in day_indices) ==
        sum(demand[t] for t in day_indices)
    )
end

# ==========================================
# 4. SOLVE
# ==========================================

optimize!(model)

# Check solve status
println("\nTermination status: ", termination_status(model))
println("Primal status: ", primal_status(model))

# ==========================================
# 5. RESULTS
# ==========================================

println("\n--- ANNUAL SYSTEM RESULTS ---")
println("Optimal Solar Size:   ", round(value(C_pv), digits=2), " kW")
println("Optimal Battery Size: ", round(value(C_batt), digits=2), " kWh")
println("Battery Power Limit:  ", round(0.25 * value(C_batt), digits=2), " kW")
println("Total Annual Cost: \$", round(objective_value(model), digits=2))

# --- COST BREAKDOWN ---
pv_capex_cost = value(C_pv) * k_pv_ann
pv_fom_cost   = value(C_pv) * fom_pv

batt_capex_cost = value(C_batt) * k_batt_ann
batt_fom_cost   = value(C_batt) * fom_batt

println("\n--- COST BREAKDOWN ---")
println("PV annualized CAPEX:      \$", round(pv_capex_cost, digits=2))
println("PV FOM:                   \$", round(pv_fom_cost, digits=2))
println("Battery annualized CAPEX: \$", round(batt_capex_cost, digits=2))
println("Battery FOM:              \$", round(batt_fom_cost, digits=2))
println("Sum of components:        \$", round(
    pv_capex_cost + pv_fom_cost +
    batt_capex_cost + batt_fom_cost, digits=2))

# ==========================================
# 6. DERIVE CHARGE / DISCHARGE FOR REPORTING
# ==========================================

p_batt_opt = value.(p_batt)

# Split signed battery power into separate series for export/plotting
battery_discharge = max.(p_batt_opt, 0.0)
battery_charge    = max.(-p_batt_opt, 0.0)

# ==========================================
# 7. ENERGY METRICS
# ==========================================

total_demand_MWh     = sum(demand[t] * dt for t in 1:T) / 1000
pv_generation_MWh    = sum(value(C_pv) * solar_30min[t] * dt for t in 1:T) / 1000
batt_charge_MWh      = sum(battery_charge[t] * dt for t in 1:T) / 1000
batt_discharge_MWh   = sum(battery_discharge[t] * dt for t in 1:T) / 1000

battery_cycles = value(C_batt) > 0 ?
    batt_discharge_MWh * 1000 / value(C_batt) : 0

println("\n--- ENERGY METRICS ---")
println("Total demand:           ", round(total_demand_MWh, digits=2), " MWh")
println("PV generation:          ", round(pv_generation_MWh, digits=2), " MWh")
println("Battery charge energy:  ", round(batt_charge_MWh, digits=2), " MWh")
println("Battery discharge:      ", round(batt_discharge_MWh, digits=2), " MWh")
println("Approx battery cycles:  ", round(battery_cycles, digits=2))

# ==========================================
# 8. QUICK CHECK: NO SIMULTANEOUS CHARGE / DISCHARGE
# ==========================================

simul = sum((battery_charge[t] > 1e-6) && (battery_discharge[t] > 1e-6) for t in 1:T)
println("\nTimesteps with simultaneous charge/discharge: ", simul)

# ==========================================
# 9. EXPORT OPTIMISED DEMAND
# ==========================================

l_actual_opt = value.(l_actual)

df_demand = DataFrame(
    timestep = 1:T,
    original_demand_kW = demand,
    optimised_demand_kW = l_actual_opt
)

CSV.write("C:/Users/naeem/Downloads/optimised_demand.csv", df_demand)

println("\nOptimised demand exported to optimised_demand.csv")

# ==========================================
# 10. EXPORT SYSTEM DISPATCH FOR MATLAB
# ==========================================

df_dispatch = DataFrame(
    timestep = 1:T,
    original_demand_kW = demand,
    optimised_demand_kW = l_actual_opt,
    pv_generation_kW = value(C_pv) .* solar_30min,
    battery_charge_kW = battery_charge,
    battery_discharge_kW = battery_discharge,
    battery_net_power_kW = p_batt_opt,
    battery_energy_kWh = value.(e_batt)
)

CSV.write("C:/Users/naeem/Downloads/system_dispatch.csv", df_dispatch)

println("Full dispatch exported to system_dispatch.csv")