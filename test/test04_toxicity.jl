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

defaultparams.spc.k_D_j
defaultparams.spc.e_z
defaultparams.spc.b_z

@testset "Toxicity with default DEB parameters" begin 
    global p = deepcopy(defaultparams)

    pmoa_idx = 6

    p.glb.t_max = 60.
    #p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = 15.

    p.spc.H_j1 = 0.2
    p.spc.H_p = 40.

    p.spc.k_D_j[pmoa_idx] = 1.
    p.spc.e_z[pmoa_idx] = 2. 
    p.spc.b_z[pmoa_idx] = 2.

    global sim = exposure(p -> ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe, 
            alg = Tsit5()
            ), 
            p,
            Matrix(hcat([0.; 1.; 2.; 4.;]...)')
            )
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :C_W_1),
        plot(:t, :E_mt, group = :C_W_1),
        plot(:t, :y_j_1_1, group = :C_W_1),
        plot(:t, :y_j_1_2, group = :C_W_1),
        plot(:t, :y_j_1_3, group = :C_W_1),
        plot(:t, :y_j_1_4, group = :C_W_1),
        plot(:t, :R, group = :C_W_1),
        plot(:t, :H, group = :C_W_1),
        xrotation = 45, 
        xlabel = "t", ylabel = ["S" "E_mt" "y_G" "y_M" "y_A" "y_R" "R" "H"], 
        leg = [:topleft false false false false], legendtitle = "C_W", 
        size = (1000,600), bottommargin = 5mm, leftmargin = 5mm
        )
    display(plt)

    #TODO: add proper test conditions

    #@test 55 <= maximum(sim.S) <= 60 # check final structural mass
    #@test ([isapprox(1, sum([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult])) for r in eachrow(sim)] |> unique)==[1] # check that exactly one life stage at a time is "true"
end

@testset "PMoA kappa" begin 
    global p = deepcopy(defaultparams)
<
    pmoa_idx = 6

    p.glb.t_max = 60.
    #p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = 15.

    p.spc.H_j1 = 0.2
    p.spc.H_p = 40.

    p.spc.k_D_j[pmoa_idx] = 1.
    p.spc.e_z[pmoa_idx] = 2. 
    p.spc.b_z[pmoa_idx] = 2.

    global sim = exposure(p -> ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe, 
            alg = Tsit5()
            ), 
            p,
            Matrix(hcat([0.; 1.; 2.; 4.;]...)')
            )
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :C_W_1),
        plot(:t, :E_mt, group = :C_W_1),
        plot(:t, :y_j_1_3, group = :C_W_1),
        plot(:t, :R, group = :C_W_1),
        plot(:t, :H, group = :C_W_1),
        xrotation = 45, 
        xlabel = "t", ylabel = ["S" "E_mt" "y_A" "R" "H"], 
        leg = [:topleft false false false false], legendtitle = "C_W", 
        size = (1000,600), bottommargin = 5mm, leftmargin = 5mm
        )
    display(plt)

    #TODO: add proper test conditions

    #@test 55 <= maximum(sim.S) <= 60 # check final structural mass
    #@test ([isapprox(1, sum([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult])) for r in eachrow(sim)] |> unique)==[1] # check that exactly one life stage at a time is "true"
end