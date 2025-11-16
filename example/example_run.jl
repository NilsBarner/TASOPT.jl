using Revise

# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
include(TASOPT.__TASOPTindices__)
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
# #  load default model
# example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# # Alternatively you can load your desired input file 
# # example_ac = read_aircraft_model("../src/IO/input.toml") # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/atr72600_input.toml")) # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/atr72600_lh2_input.toml")) # MODIFY <path> appropriately
example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/a220100_input.toml")) # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/a220100_lh2_input_org.toml")) # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml")) # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(pwd(), "../nils/point_loads/default_regional_cryo.toml")) # MODIFY <path> appropriately
# example_ac = read_aircraft_model(joinpath(pwd(), "nils", "point_loads", "atr72600_lh2_input.toml"))

# 3) Size aircraft
time_size_aircraft = @elapsed size_aircraft!(example_ac; iter=100) # second argument added by NILS
println("Time to size aircraft = $time_size_aircraft s")

# 4) Visualize outputs
# Output resulting geometry of aircraft
summary(example_ac)
# Or individually look at certain aspects:
# Show weight buildup of the aircraft:
# TASOPT.weight_buildup(example_ac) 
# # Show aerodynamics:
# TASOPT.aero(example_ac)
# # Geometry:
# TASOPT.geometry(example_ac)

# 5) Plot figures
using Plots
# p = TASOPT.stickfig(example_ac)  # line commented by Nils
p = TASOPT.stickfig(example_ac, show_grid = false, annotate_text = false, show_seats = false, fill_tank = false, annotate_sketch = false, annotate_length = false, annotate_group = false, airframe_colour = :black, engine_colour = :black)  # line added by Nils
# savefig(p, "Example.png")
savefig(p, "Example.svg")
display(p)