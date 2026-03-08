using JuMP, HiGHS, CSV, DataFrames, Interpolations, Statistics, Dates

# ==========================================
# 1. DATA LOADING & ALIGNMENT
# ==========================================

df_dc = CSV.read("D:/Masters/Imperial/Modules/Power System Planning/DATA/datacentre_nsw_merged5.csv", DataFrame)
df_grid_load = CSV.read("D:/Masters/Imperial/Modules/Power System Planning/DATA/NSW_2019_combined.csv", DataFrame)

df_solar = CSV.read("D:/Masters/Imperial/Modules/Power System Planning/DATA/irradiance_AU_2019_extracted.csv", DataFrame)
df_wind  = CSV.read("D:/Masters/Imperial/Modules/Power System Planning/DATA/ninja_wind_-37.9475_143.7511_uncorrected.csv", DataFrame, header=4)

# ---- Solar ----
raw_irradiance = df_solar.AU
solar_pu_hourly = min.(raw_irradiance ./ 1000.0, 1.0)

# ---- Wind ----
wind_pu_hourly = df_wind.electricity

# Shift 20 intervals (10 hours) to convert UTC to AEST
solar_pu_hourly = circshift(solar_pu_hourly, 20)
wind_pu_hourly = circshift(wind_pu_hourly, 20)

# ---- Time alignment ----
T_intervals = nrow(df_grid_load)
T_hours = length(wind_pu_hourly)
dt = 0.5

h_steps = 1:T_hours
i_steps = range(1, stop=T_hours, length=T_intervals)

itp_solar  = linear_interpolation(h_steps, solar_pu_hourly, extrapolation_bc=Flat())
itp_wind  = linear_interpolation(h_steps, wind_pu_hourly, extrapolation_bc=Flat())

pv_pu_30   = [itp_solar(t)  for t in i_steps]
wind_pu_30 = [itp_wind(t)  for t in i_steps]

aus_load_30 = df_grid_load.TOTALDEMAND
dc_load_30  = df_dc.DATACENTRE_DEMAND
total_load  = aus_load_30 .+ dc_load_30

# ==========================================
# 2. ECONOMIC PARAMETERS
# ==========================================

r = 0.06
crf(i,n) = (i*(1+i)^n)/((1+i)^n - 1)

# Lifetimes
life_ren = 30
life_batt = 30
life_gas = 30  

# CAPEX ($/kW or $/kWh)
capex_solar_kw = 1491.64
capex_wind_kw  = 1569.099
capex_gas_kw   = 1727.7    
capex_batt_kwh = 1711.4

# FOM ($/kW-year or $/kWh-year)
fom_solar_kw = 21.001
fom_wind_kw  = 31.251
fom_batt_kwh = 38.78
fom_gas_kw   = 38.7        

# VOM - ($/MWh)
vom_gas = 52.67            
                           
# Annualized cost ($/MW-year or $/MWh-year)
k_solar_ann = (capex_solar_kw*1000)*crf(r,life_ren) + (fom_solar_kw*1000)
k_wind_ann  = (capex_wind_kw*1000)*crf(r,life_ren) + (fom_wind_kw*1000)
k_gas_ann   = (capex_gas_kw*1000)*crf(r,life_gas)  + (fom_gas_kw*1000)  
k_batt_ann  = (capex_batt_kwh*1000)*crf(r,life_batt) + (fom_batt_kwh*1000)

# ==========================================
# 3. OPTIMIZATION MODEL 
# ==========================================

model = Model(HiGHS.Optimizer)
set_silent(model)

# Solver settings - INCREASE TIME LIMIT
set_time_limit_sec(model, 1800.0)  # 30 minutes instead of 10
set_optimizer_attribute(model, "mip_rel_gap", 0.005)  # Tighter tolerance

# Investment variables
@variable(model, Cap_Solar >= 0)
@variable(model, Cap_Wind  >= 0)
@variable(model, Cap_Batt  >= 0)  
@variable(model, Cap_Gas   >= 0)  

# Operational variables
@variable(model, p_solar[1:T_intervals] >= 0)
@variable(model, p_wind[1:T_intervals]  >= 0)
@variable(model, p_gas[1:T_intervals]   >= 0)  
@variable(model, p_ch[1:T_intervals]    >= 0)
@variable(model, p_dis[1:T_intervals]   >= 0)
@variable(model, e_batt[1:T_intervals]  >= 0)

# Energy Not Served
@variable(model, ens[t=1:T_intervals] >= 0)

# DC flexibility variables
@variable(model, p_dc[1:T_intervals] >= 0)

# ==========================================
# *** RELIABILITY CONSTRAINTS (LOLP) ***
# ==========================================

@variable(model, loss_of_load[t=1:T_intervals], Bin)

# Big-M constraint
M = maximum(total_load) * 2

big_m_constraints = [@constraint(model, ens[t] <= M * loss_of_load[t]) for t in 1:T_intervals]

# LOLE Constraint
MAX_LOLE_INTERVALS = 20  # 10 hours/year

LOLE_LIMIT = @constraint(model, sum(loss_of_load[t] for t in 1:T_intervals) <= MAX_LOLE_INTERVALS)

# EUE Constraint - THIS IS CRITICAL!
total_annual_demand = sum(total_load) * dt
MAX_EUE_MWH = total_annual_demand * 0.0001  # 0.01%
EUE_LIMIT = @constraint(model, sum(ens[t] * dt for t in 1:T_intervals) <= MAX_EUE_MWH)

# Value of Lost Load
voll = 2000000

@objective(model, Min,
    Cap_Solar*k_solar_ann +
    Cap_Wind*k_wind_ann +
    Cap_Batt*k_batt_ann +
    Cap_Gas*k_gas_ann +  
    sum(p_gas[t]*vom_gas*dt for t in 1:T_intervals) +
    sum(ens[t] * voll * dt for t in 1:T_intervals)
)

# Energy balance
@constraint(model, balance[t=1:T_intervals],
    p_solar[t] + p_wind[t] + p_gas[t] + p_dis[t] + ens[t] == aus_load_30[t] + p_dc[t] + p_ch[t]
)

# DC flexibility parameters
flex_limit = 0.10
T_day = 48

for t in 1:T_intervals
    @constraint(model, p_dc[t] >= dc_load_30[t] * (1 - flex_limit))
    @constraint(model, p_dc[t] <= dc_load_30[t] * (1 + flex_limit))
end

# Energy Neutrality
for d in 1:Int(T_intervals/T_day)
    day_indices = ( (d-1)*T_day + 1 ) : (d*T_day)
    @constraint(model, sum(p_dc[t] for t in day_indices) == sum(dc_load_30[t] for t in day_indices))
end

# *** FIXED GAS CONSTRAINT ***
# Gas capacity must be ≤ 2% of total generation capacity
@constraint(model, gas_limit, 
    Cap_Gas * 0.98 <= 0.02 * (Cap_Solar + Cap_Wind + Cap_Batt*0.25))

# Battery efficiency
eff = 0.95

# Renewable generation limits
@constraint(model, solar_limit[t=1:T_intervals], p_solar[t] <= Cap_Solar * pv_pu_30[t])
@constraint(model, wind_limit[t=1:T_intervals], p_wind[t]  <= Cap_Wind  * wind_pu_30[t])

# Battery and gas constraints
for t in 1:T_intervals
    # Gas generation limit
    @constraint(model, p_gas[t] <= Cap_Gas)
    
    # Battery power limits (4-hour system)
    @constraint(model, p_ch[t]  <= Cap_Batt * 0.25)
    @constraint(model, p_dis[t] <= Cap_Batt * 0.25)

    # Battery energy dynamics
    if t == 1
        @constraint(model, e_batt[t] == e_batt[T_intervals]
            + p_ch[t]*dt*eff - p_dis[t]*dt/eff)
    else
        @constraint(model, e_batt[t] == e_batt[t-1]
            + p_ch[t]*dt*eff - p_dis[t]*dt/eff)
    end

    # Battery energy capacity limit
    @constraint(model, e_batt[t] <= Cap_Batt)
end

println("\n✓ Model built successfully")
println("  Variables: ", num_variables(model))
println("  Constraints: ", num_constraints(model; count_variable_in_set_constraints=false))

# ==========================================
# 4. TWO-PASS SOLUTION
# ==========================================

println("\n" * "="^70)
println("=== PASS 1: SOLVING MILP FOR OPTIMAL CAPACITY ===")
println("="^70)

optimize!(model)

# Check if solution is optimal
if termination_status(model) != MOI.OPTIMAL
    println("WARNING: Pass 1 did not reach optimality!")
    println("Status: ", termination_status(model))
end

# Capture optimal capacities
opt_cap_solar = value(Cap_Solar)
opt_cap_wind  = value(Cap_Wind)
opt_cap_gas   = value(Cap_Gas)  
opt_cap_batt  = value(Cap_Batt)
opt_cap_batt_power = opt_cap_batt * 0.25  # Convert to power (MW)

# Calculate gas percentage
total_capacity = opt_cap_solar + opt_cap_wind + opt_cap_batt_power + opt_cap_gas
gas_percentage = 100 * opt_cap_gas / total_capacity

println("\nPass 1 Complete - Optimal Capacities:")
println("  Solar:          ", round(opt_cap_solar, digits=2), " MW")
println("  Wind:           ", round(opt_cap_wind, digits=2), " MW")
println("  Gas:            ", round(opt_cap_gas, digits=2), " MW")
println("  Battery Energy: ", round(opt_cap_batt, digits=2), " MWh")
println("  Battery Power:  ", round(opt_cap_batt_power, digits=2), " MW")
println("  Gas Percentage: ", round(gas_percentage, digits=2), " %")

# Check Pass 1 reliability
pass1_ens = sum(value.(ens)) * dt
pass1_eue_pct = 100 * pass1_ens / total_annual_demand
println("\n  Pass 1 EUE: ", round(pass1_ens, digits=2), " MWh (", round(pass1_eue_pct, digits=4), "%)")

if pass1_eue_pct > 0.01
    println(" WARNING: EUE exceeds 0.01% in Pass 1! Solver may need more time.")
end

# --- PASS 2: Re-solve as LP for shadow prices ---
println("\n" * "="^70)
println("=== PASS 2: RE-SOLVING AS LP FOR SHADOW PRICES ===")
println("="^70)

# Fix capacities
fix(Cap_Solar, opt_cap_solar; force=true)
fix(Cap_Wind,  opt_cap_wind;  force=true)
fix(Cap_Gas,   opt_cap_gas;   force=true)
fix(Cap_Batt,  opt_cap_batt;  force=true)

# Relax binary variables
unset_binary.(loss_of_load)
set_lower_bound.(loss_of_load, 0.0)
set_upper_bound.(loss_of_load, 1.0)

# Only remove LOLE constraint (which uses binaries)
delete(model, LOLE_LIMIT)

# Remove Big-M constraints (they use binaries)
for t in 1:T_intervals
    delete(model, big_m_constraints[t])
end

# Remove gas limit constraint (capacity already fixed)
delete(model, gas_limit)

optimize!(model)

println(" Pass 2 Complete")

status = termination_status(model)
solve_time_sec = solve_time(model)

println("\nTotal solve time: ", round(solve_time_sec, digits=2), " seconds")
println("Status: ", status)

# ==========================================
# 5. EXTRACT SHADOW PRICES AND COSTS
# ==========================================

println("\n" * "="^70)
println("ELECTRICITY PRICING ANALYSIS")
println("="^70)

# Extract shadow prices from energy balance constraint
shadow_prices_marginal = dual.(balance) ./ dt

println("\n--- Shadow Price Statistics (Marginal Cost) ---")
println("Min:    \$", round(minimum(shadow_prices_marginal), digits=2), "/MWh")
println("Max:    \$", round(maximum(shadow_prices_marginal), digits=2), "/MWh")
println("Mean:   \$", round(mean(shadow_prices_marginal), digits=2), "/MWh")
println("Median: \$", round(median(shadow_prices_marginal), digits=2), "/MWh")

# Count price spikes
price_spikes = sum(shadow_prices_marginal .> 500)
println("Intervals with price > \$500/MWh: ", price_spikes, " (", round(100*price_spikes/T_intervals, digits=2), "%)")

if maximum(shadow_prices_marginal) > 1000
    println("\n WARNING: Extreme price spikes detected!")
    println("  This indicates capacity scarcity and ENS events")
end

println("\n" * "="^70)
println("Full price = Marginal price + Fixed Cost Recovery")
println()

# Calculate total fixed costs
fixed_cost_solar = opt_cap_solar * k_solar_ann
fixed_cost_wind = opt_cap_wind * k_wind_ann
fixed_cost_gas = opt_cap_gas * k_gas_ann
fixed_cost_batt = opt_cap_batt * k_batt_ann

total_fixed_cost = fixed_cost_solar + fixed_cost_wind + fixed_cost_gas + fixed_cost_batt

println("--- Fixed Costs (Annual) ---")
println("Solar:          \$", round(fixed_cost_solar/1e6, digits=2), " million/year")
println("Wind:           \$", round(fixed_cost_wind/1e6, digits=2), " million/year")
println("Gas:            \$", round(fixed_cost_gas/1e6, digits=2), " million/year")
println("Battery:        \$", round(fixed_cost_batt/1e6, digits=2), " million/year")
println("Total Fixed:    \$", round(total_fixed_cost/1e6, digits=2), " million/year")

# Total energy delivered
total_energy_delivered = sum(total_load) * dt

# Fixed Cost Recovery
fixed_cost_recovery = total_fixed_cost / total_energy_delivered

println("\n--- Fixed Cost Recovery ---")
println("Total energy delivered: ", round(total_energy_delivered/1e6, digits=2), " million MWh/year")
println("Fixed cost per MWh:     \$", round(fixed_cost_recovery, digits=2), "/MWh")

# Full electricity price
retail_price = shadow_prices_marginal .+ fixed_cost_recovery

println("\n--- Full Electricity Price (Retail) ---")
println("Retail price = Marginal price + Fixed cost recovery")
println("Min:    \$", round(minimum(retail_price), digits=2), "/MWh")
println("Max:    \$", round(maximum(retail_price), digits=2), "/MWh")
println("Mean:   \$", round(mean(retail_price), digits=2), "/MWh")

println("\n" * "="^70)

# DC costs
dc_energy_consumed = sum(value.(p_dc)) * dt
dc_bill_marginal = sum(value.(p_dc) .* shadow_prices_marginal .* dt)

println("\n--- Marginal Pricing (Short-term cost) ---")
println("DC energy consumed: ", round(dc_energy_consumed/1e6, digits=3), " million MWh/year")
println("DC bill (marginal): \$", round(dc_bill_marginal/1e6, digits=2), " million/year")

dc_bill_retail = sum(value.(p_dc) .* retail_price .* dt)

println("\n--- Retail Pricing (Full cost recovery) ---")
println("DC bill (retail):   \$", round(dc_bill_retail/1e6, digits=2), " million/year")

# ==========================================
# 6. RELIABILITY METRICS
# ==========================================

ens_vals = value.(ens)
lol_vals = value.(loss_of_load)

LOLE_intervals = sum(lol_vals .> 0.5)
LOLE_hours = LOLE_intervals * dt
EUE_MWh = sum(ens_vals .* dt)
EUE_percentage = 100 * EUE_MWh / total_annual_demand
LOLP = LOLE_intervals / T_intervals

println("\n" * "="^70)
println("RELIABILITY METRICS")
println("="^70)
println("LOLE: ", round(LOLE_hours, digits=2), " hours/year")
println("EUE:  ", round(EUE_MWh, digits=4), " MWh/year (", round(EUE_percentage, digits=6), "%)")
println("LOLP: ", round(LOLP*100, digits=4), "%")

# Verification
println("\n--- Reliability Verification ---")
if EUE_percentage <= 0.01
    println(" EUE constraint MET (", round(EUE_percentage, digits=6), "% ≤ 0.01%)")
else
    println(" EUE constraint VIOLATED (", round(EUE_percentage, digits=6), "% > 0.01%)")
    println("  This should not happen - check solver settings!")
end

# Check battery operation
p_ch_vals = value.(p_ch)
p_dis_vals = value.(p_dis)

simultaneous = sum((p_ch_vals .> 0.01) .& (p_dis_vals .> 0.01))

println("\n--- BATTERY OPERATION CHECK ---")
println("Intervals with simultaneous charge/discharge: ", simultaneous)

if simultaneous > 0
    println("Battery charging and discharging simultaneously in ", simultaneous, " intervals")
    max_simultaneous = maximum(min.(p_ch_vals, p_dis_vals))
    println("  Max simultaneous power: ", round(max_simultaneous, digits=2), " MW")
else
    println(" No simultaneous charge/discharge detected")
end

# ==========================================
# 7. GENERATION STATISTICS
# ==========================================

gen_solar = sum(value.(p_solar)) * dt
gen_wind = sum(value.(p_wind)) * dt
gen_gas = sum(value.(p_gas)) * dt
gen_total = gen_solar + gen_wind + gen_gas

println("\n" * "="^70)
println("GENERATION STATISTICS")
println("="^70)
println("Solar:  ", round(gen_solar/1e6, digits=3), " million MWh (", round(100*gen_solar/gen_total, digits=1), "%)")
println("Wind:   ", round(gen_wind/1e6, digits=3), " million MWh (", round(100*gen_wind/gen_total, digits=1), "%)")
println("Gas:    ", round(gen_gas/1e6, digits=3), " million MWh (", round(100*gen_gas/gen_total, digits=1), "%)")
println("Total:  ", round(gen_total/1e6, digits=3), " million MWh")

# *** CRITICAL: CALCULATE BATTERY STATISTICS HERE ***
batt_charge_energy = sum(p_ch_vals) * dt
batt_discharge_energy = sum(p_dis_vals) * dt
batt_efficiency_realized = batt_discharge_energy / max(batt_charge_energy, 1.0) * 100
batt_cycles = batt_discharge_energy / max(opt_cap_batt, 1.0)

println("\n--- Battery Performance ---")
println("Charged:        ", round(batt_charge_energy/1e6, digits=3), " million MWh")
println("Discharged:     ", round(batt_discharge_energy/1e6, digits=3), " million MWh")
println("Efficiency:     ", round(batt_efficiency_realized, digits=1), "%")
println("Annual cycles:  ", round(batt_cycles, digits=1), " cycles/year")

# LCOE calculations
lcoe_solar = (opt_cap_solar * k_solar_ann) / max(gen_solar, 1.0)
lcoe_wind = (opt_cap_wind * k_wind_ann) / max(gen_wind, 1.0)
lcoe_gas = ((opt_cap_gas * k_gas_ann) + (gen_gas * vom_gas)) / max(gen_gas, 1.0)

total_system_cost = total_fixed_cost + (gen_gas * vom_gas)
system_lcoe = total_system_cost / total_energy_delivered

println("\n--- LCOE (Levelized Cost of Energy) ---")
println("Solar:  \$", round(lcoe_solar, digits=2), "/MWh")
println("Wind:   \$", round(lcoe_wind, digits=2), "/MWh")
println("Gas:    \$", round(lcoe_gas, digits=2), "/MWh")
println("System: \$", round(system_lcoe, digits=2), "/MWh")

# ==========================================
# 8. EXPORT RESULTS
# ==========================================

println("\n" * "="^70)
println("EXPORTING RESULTS")
println("="^70)

# Battery SoC
soc_percent = if opt_cap_batt > 0 
    [value(e_batt[t]) / opt_cap_batt * 100 for t in 1:T_intervals]
else
    zeros(T_intervals)
end

# Timestamps
start_datetime = DateTime(2019, 1, 1, 11, 0, 0)
timestamps = [start_datetime + Minute(30*(i-1)) for i in 1:T_intervals]

# Create results DataFrame
results_df = DataFrame(
    Timestamp = timestamps,
    Shadow_Price_Marginal_MWh = shadow_prices_marginal,  
    Retail_Price_MWh = retail_price,                      
    Fixed_Cost_Recovery_MWh = fill(fixed_cost_recovery, T_intervals),
    DC_Demand_Actual_MW = value.(p_dc),
    DC_Demand_Baseline_MW = dc_load_30,
    Grid_Demand_MW = aus_load_30,
    Total_Demand_MW = total_load,
    Solar_Gen_MW = value.(p_solar),
    Wind_Gen_MW = value.(p_wind),
    Gas_Gen_MW = value.(p_gas),
    Batt_Discharge_MW = value.(p_dis),
    Batt_Charge_MW = value.(p_ch),
    Battery_SoC_Percent = soc_percent,
    ENS_MW = ens_vals,
    Loss_of_Load_Flag = lol_vals
)

CSV.write("D:/Masters/Imperial/Modules/Power System Planning/RESULTS/Grid_Level_Optimization_Results.csv", results_df)
println(" Main results exported to: Grid_Level_Optimization_Results.csv")

# Summary metrics
summary_df = DataFrame(
    Metric = [
        "Solar Capacity (MW)",
        "Wind Capacity (MW)",
        "Gas Capacity (MW)",
        "Battery Power (MW)",
        "Battery Energy (MWh)",
        "Battery Duration (hours)",
        "Gas % of Total Capacity",
        "LOLE (hours/year)",
        "EUE (MWh/year)",
        "EUE (% of demand)",
        "LOLP (%)",
        "Max Shadow Price (USD/MWh)",
        "Mean Shadow Price (USD/MWh)",
        "Fixed Cost Recovery (USD/MWh)",
        "Mean Retail Price (USD/MWh)",
        "System LCOE (USD/MWh)",
        "DC Bill Marginal (USD million/year)",
        "DC Bill Retail (USD million/year)",
        "Total System Cost (USD million/year)",
        "Solar Generation (million MWh)",
        "Wind Generation (million MWh)",
        "Gas Generation (million MWh)",
        "Battery Charge (million MWh)",
        "Battery Discharge (million MWh)",
        "Battery Efficiency (%)",
        "Battery Cycles (per year)",
        "Gas VOM (USD/MWh)"
    ],
    Value = [
        opt_cap_solar,
        opt_cap_wind,
        opt_cap_gas,
        opt_cap_batt_power,
        opt_cap_batt,
        4.0,
        gas_percentage,
        LOLE_hours,
        EUE_MWh,
        EUE_percentage,
        LOLP*100,
        maximum(shadow_prices_marginal),
        mean(shadow_prices_marginal),
        fixed_cost_recovery,
        mean(retail_price),
        system_lcoe,
        dc_bill_marginal/1e6,
        dc_bill_retail/1e6,
        total_system_cost/1e6,
        gen_solar/1e6,
        gen_wind/1e6,
        gen_gas/1e6,
        batt_charge_energy/1e6,
        batt_discharge_energy/1e6,
        batt_efficiency_realized,
        batt_cycles,
        vom_gas
    ]
)

CSV.write("D:/Masters/Imperial/Modules/Power System Planning/RESULTS/Grid_Level_Optimization_Summary.csv", summary_df)
println(" Summary exported to: Grid_Level_Optimization_Summary.csv")

println("\n" * "="^70)
println("OPTIMIZATION COMPLETE!")
println("="^70)
println("\nKey Results:")
println("  Solar:          ", round(opt_cap_solar, digits=2), " MW")
println("  Wind:           ", round(opt_cap_wind, digits=2), " MW")
println("  Gas:            ", round(opt_cap_gas, digits=2), " MW (", round(gas_percentage, digits=2), "%)")
println("  Battery Energy: ", round(opt_cap_batt, digits=2), " MWh")
println("  Battery Power:  ", round(opt_cap_batt_power, digits=2), " MW")
println("  EUE:            ", round(EUE_MWh, digits=2), " MWh (", round(EUE_percentage, digits=4), "%)")
println("  Max price:      \$", round(maximum(shadow_prices_marginal), digits=2), "/MWh")
println("  Mean price:     \$", round(mean(shadow_prices_marginal), digits=2), "/MWh")
println("  System LCOE:    \$", round(system_lcoe, digits=2), "/MWh")
println("  DC bill (retail): \$", round(dc_bill_retail/1e6, digits=2), " million/year")
println("="^70)
