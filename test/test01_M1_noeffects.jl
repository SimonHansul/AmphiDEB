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
    p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = 15.

    p.spc.H_j1 = 1.

    @time global sim = ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe, 
            saveat = 1/240,
            alg = Tsit5()
            );

    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:dI] = vcat(0, diff(sim.I))
    sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = plot_statevars(
        sim, 
        [:S, :H, :E_mt_rel, :R, :X_emb, :J, :dI_rel, :W_tot, :f_X], 
        xrotation = 45
        )
    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)

    display(plt)

    # check final structural mass
    @test 6 <= maximum(sim.S) <= 8 
    # check sum of life stage indicators - can deviate a little but not much from 1
    #@test ([isapprox(1, sum([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult]), atol = 1e-3) for r in eachrow(sim)] |> unique)==[1] # check that exactly one life stage at a time is "true"
end

using BenchmarkTools
@benchmark ODE_simulator(
    p, 
    returntype = EcotoxSystems.dataframe, 
    alg = nothing
    )

