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

@testset "Pathogen default parameters: infection without effects" begin 
    global p = deepcopy(defaultparams)

    p.glb.t_max = 60
    p.glb.pathogen_inoculation_time = 20.
    p.glb.pathogen_inoculation_dose = 1e3

    p.glb.dX_in = [15., 15.]

    p.spc.Z = Dirac(1.)

    global sim = @replicates ODE_simulator(
            p, 
            returntype = EcotoxSystems.dataframe
            ) 10
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :replicate),
        plot(:t, :E_mt, group = :replicate),
        plot(:t, :P_Z, group = :replicate),
        xrotation = 45
        )
    
    @df sim lineplot!(:t, :P_Z, subplot = 3, lw = 5, fillalpha = .2)

    # check for mean and variability in pathogen loads

    let df = @subset(sim, :t .== maximum(:t))

        P_Z_mean = mean(df.P_Z)
        P_Z_cv = std(df.P_Z) / P_Z_mean

        @test 500 < P_Z_mean < 1500
        @test 0.5 < P_Z_cv < 2
    end
    
    display(plt)
end


@testset "Pathogen default parameters: infection with effects" begin 

    pmoa_P_idx = 1

    global p = deepcopy(defaultparams)

    p.glb.t_max = 60
    p.glb.pathogen_inoculation_time = 20.
    p.glb.pathogen_inoculation_dose = 1e3

    p.glb.dX_in = [15., 15.]

    p.spc.Z = Dirac(1.)
    p.spc.E_P[pmoa_P_idx] = 1. 
    p.spc.B_P[pmoa_P_idx] = 2.

    @time global sim = @replicates ODE_simulator(p) 10
        
    sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
    sim[!,:W_tot] = sim.S .+ sim.E_mt 
    
    plt = @df sim plot(
        plot(:t, :S, group = :replicate, ylabel = "S"),
        plot(:t, :E_mt, group = :replicate, ylabel = "E_mt"),
        plot(:t, :P_S, group = :replicate, ylabel = "P_S"),
        plot(:t, :P_Z, group = :replicate, ylabel = "P_Z"),
        xrotation = 45, xlabel = "t"
        )

    @df sim lineplot!(plt, :t, :P_S, lw = 5, fillalpha = .2, subplot = 3)
    @df sim lineplot!(plt, :t, :P_Z, lw = 5, fillalpha = .2, subplot = 4)
    
    # check for mean and variability in pathogen loads

    let df = @subset(sim, :t .== maximum(:t))

        P_Z_mean = mean(df.P_Z)
        P_Z_cv = std(df.P_Z) / P_Z_mean

        @test 500 < P_Z_mean < 1500
        @test 0.25 < P_Z_cv < 2
    end
    
    display(plt)
end


