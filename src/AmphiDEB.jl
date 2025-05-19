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
include("derivatives.jl") # derivatives of the default model
include("statevars.jl") # setting up state variables
include("individual_rules.jl") # individual rule-based component
include("global_rules.jl") # global rule-based component
include("simulators.jl") # functions to run run simulations

include("traits.jl") # functions to infer traits from parameters or simulation output (e.g. maximum size, age at birth, etc.)
include("utils.jl") # various auxiliary functions

# to precompile the model, we simulate the default parameters
@compile_workload begin
    p = deepcopy(defaultparams)

    sim = @replicates ODE_simulator(p) 10
    
    p.glb.t_max = 365.
    p.glb.dX_in = [500., 500.]
    p.spc.tau_R = 30.
    
    sim = IBM_simulator(p)
end

end # module AmphiDEB
