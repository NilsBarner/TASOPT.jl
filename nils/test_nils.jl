# Use GitHub dev installation
# import Pkg
# Pkg.develop(path="C:\\Users\\nmb48\\Documents\\GitHub\\TASOPT.jl")

# Revert back to simple installation
# Pkg.free("TASOPT")

using TASOPT
example_ac = load_default_model()
size_aircraft!(example_ac)