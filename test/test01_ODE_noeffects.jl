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


@testset "Default parameters" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365*2
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = 20.
    p.spc.H_p = 50.

    S_max_anl = AmphiDEB.calc_S_max_juv(p.spc)

    @time global sim = AmphiDEB.ODE_simulator(
            p, 
            saveat = 1,
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
    @test 0.8*S_max_anl <= maximum(sim.S) <= 1.2*S_max_anl 

    # check that all life stage indicators max out close to 1
    @test 0.99 < maximum(sim.embryo) < 1.01
    @test 0.99 < maximum(sim.larva) < 1.01
    @test 0.99 < maximum(sim.metamorph) < 1.01
    @test 0.99 < maximum(sim.juvenile) < 1.01
    @test 0.99 < maximum(sim.adult) < 1.01

    # check that the sum of life stage indicators is always approximately 1

    sum_indicators = @. sim.embryo + sim.larva + sim.metamorph + sim.juvenile + sim.adult

    @test unique(isapprox.(1, sum_indicators, atol = 1e-3)) == [true]
end



@testset "Randomized parameters" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365*2
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = 20.


    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.k_M_emb = truncated(Normal(0.11, 0.011), 0, Inf)
    p.spc.eta_AR = truncated(Normal(0.95, 0.095), 0, 1)

    p.spc.H_p = 50.

    S_max_anl = AmphiDEB.calc_S_max_juv(p.spc)

    @time global sim = @replicates AmphiDEB.ODE_simulator(
            p, 
            saveat = 1,
            alg = Tsit5()
            ) 10

    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:dI] = vcat(0, diff(sim.I))
    sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        lineplot(:t, :S), 
        lineplot(:t, :H), 
        lineplot(:t, :E_mt_rel), 
        lineplot(:t, :R),
        xrotation = 45
        )
    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)

    display(plt)
    
    # check that all life stage indicators max out close to 1
    @test 0.99 < maximum(sim.embryo) < 1.01
    @test 0.99 < maximum(sim.larva) < 1.01
    @test 0.99 < maximum(sim.metamorph) < 1.01
    @test 0.99 < maximum(sim.juvenile) < 1.01
    @test 0.99 < maximum(sim.adult) < 1.01

    # check that the sum of life stage indicators is always approximately 1

    sum_indicators = @. sim.embryo + sim.larva + sim.metamorph + sim.juvenile + sim.adult

    @test unique(isapprox.(1, sum_indicators, atol = 1e-3)) == [true]
end
