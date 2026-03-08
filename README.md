# Datacentre-Energy-Optimisation-and-PV-Battery-Microgrid-Model-using-Julia
In this project, Julia is used to model the cost of energy for a Datacentre in New South Wales, Australia. Their energy cost when buying from the grid is modelled and compared to the potential cost of independently investing in a solar and BESS system to see if it is a wise investment.

This Github repo contaiins all the datasets needed to run both models as well as the model codes used. You can download the datasets and replace the paths in the codes available.

We found it to be substantially cheaper (5.5 billion dollars cheaper)for the datacentre to remain connected to the grid.

You can read about our method and assumptions in the "Data Centre Optimisation Report"

In the "Datasets" branch,

-30-minute Data Centre demand data = datacentre_nsw_merged5.csv

-hourly solar data = ninja_pv_-31.3665_144.8481_uncorrected.csv

In the "Codes" section,

-The model for the data independent of the gris is "datacentre_optimization_nogrid.jl"

