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

# Amphibian model without chemical or pathogen effects 
Amphibian_noeffects!(du, u, p, t) = begin
    DEBODE_global!(du, u, p, t)
    u.ind.y_j .= 1 # making sure that all effects are turned off
    u.ind.y_jP .= 1
    Amphibian_DEB_M1!(du, u, p, t)
end

# implementation works as long as k_J is not too high...
# next: trying a gradual change to metamorph

defaultparams.spc.dI_max_juv = 1
AmphiDEB.calc_S_max_juv(defaultparams.spc)

@testset "Default parameters" begin 
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365*2
    p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = 15.

    @time global sim = ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe, 
            alg = Tsit5()
            );

    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:dI] = vcat(0, diff(sim.I))
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = plot_statevars(
        sim, 
        [:S, :H, :E_mt_rel, :R, :X_emb, :J, :R, :dI, :W_tot, :f_X], 
        xrotation = 45
        )
    hline!([p.spc.H_j1], subplot=2, color = :gray, linestyle = :dash)

    display(plt)

    @test 55 <= maximum(sim.S) <= 60 # check final structural mass
    @test ([sum([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult])==1 for r in eachrow(sim)] |> unique)==[1] # check that exactly one life stage at a time is "true"
end


