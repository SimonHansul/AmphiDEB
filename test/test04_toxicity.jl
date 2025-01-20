using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta
using StatsBase

using Revise

@time import AmphiDEB: defaultparams, ODE_simulator, Amphibian_DEB!, AmphiDEB_ODE!
using AmphiDEB
using EcotoxSystems
import EcotoxSystems: DEBODE_global!

import EcotoxSystems: sig
import EcotoxSystems: constrmvec

defaultparams.spc.dI_max_juv = 1
AmphiDEB.calc_S_max_juv(defaultparams.spc)

@testset "Default parameters" begin 
    global p = deepcopy(defaultparams)

    p.glb.t_max = 60
    #p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = 15.

    p.spc.k_D_z[3] = 1.
    p.spc.e_z[3] = 1. 
    p.spc.b_z[3] = .1

    sim = exposure(p -> ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe, 
            alg = Tsit5()
            ), 
            p,
            Matrix(hcat([0.; 0.5; 1.5; 2; 4.]...)')
            )
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :C_W_1),
        plot(:t, :E_mt, group = :C_W_1),
        xrotation = 45, 
        xlabel = "t", ylabel = ["S" "E_mt"], leg = [:topleft false], legendtitle = "C_W"
        )
    display(plt)

    #TODO: add proper test conditions

    #@test 55 <= maximum(sim.S) <= 60 # check final structural mass
    #@test ([isapprox(1, sum([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult])) for r in eachrow(sim)] |> unique)==[1] # check that exactly one life stage at a time is "true"
end