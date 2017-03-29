include("Include.jl")

# load the data dictionary -
data_dictionary = DataDictionary(0,0,0)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
