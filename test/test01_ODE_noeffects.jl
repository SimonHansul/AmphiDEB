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

    @test 0.8*S_max_anl <= maximum(sim.S) <= 1.2*S_max_anl # check final structural mass
end
