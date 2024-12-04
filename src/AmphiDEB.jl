module AmphiDEB

using PrecompileTools
using EcotoxSystems
import EcotoxSystems: sig, clipneg

using Parameters
using Distributions
using DataStructures
import DataFrames: AbstractDataFrame
using ComponentArrays, StaticArrays
using OrdinaryDiffEq
using StatsBase

include("default_params.jl") # default parameter sets
include("derivatives.jl") # definition of model derivatives
include("statevars.jl") # setting up state variables
include("simulators.jl") # running simulations

include("traits.jl") # functions to infer traits from parameters or simulation output (e.g. maximum size, age at birth, etc.)
include("utils.jl") # various auxiliary functions
export plot_statevars

#@compile_workload begin
#    yhat = ODE_simulator(defaultparams)
#end


end # module AmphiDEB
