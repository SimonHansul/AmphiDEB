using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta
using StatsBase

using EcotoxSystems
import EcotoxSystems: constrmvec

using Revise

using AmphiDEB
@time import AmphiDEB: defaultparams, ODE_simulator, Amphibian_DEB!, AmphiDEB_ODE!


@testset "Toxicity with default parameters" begin 
    
    for pmoa_idx in 1:6
    
        global p = deepcopy(defaultparams)
    
        p.glb.t_max = 60.
        #p.glb.pathogen_inoculation_time = Inf

        p.glb.dX_in = [15. 15.]

        p.spc.H_j1 = 0.2
        p.spc.H_p = 40.

        p.spc.KD .= 0.
        p.spc.KD[pmoa_idx] = 1.
        p.spc.E[pmoa_idx] = 2. 
        p.spc.B[pmoa_idx] = 2.

        # check that the syntax for parameter assignment works 
        @test sum(p.spc.KD .> 0) == 1

        @time global sim = EcotoxSystems.exposure(
            ODE_simulator, 
            p,
            Matrix(hcat([0.; 1.; 2.;]...)')
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
            plot(:t, :y_j_1_5, group = :C_W_1),
            plot(:t, :y_j_1_6, group = :C_W_1),
            plot(:t, :R, group = :C_W_1),
            plot(:t, :H, group = :C_W_1),
            xrotation = 45, 
            xlabel = "t", ylabel = ["S" "E_mt" "y_G" "y_M" "y_A" "y_R" "y_H" "y_κ" "R" "H"], 
            leg = [:topleft false false false false], legendtitle = "C_W", 
            size = (1000,600), bottommargin = 5mm, leftmargin = 5mm
            )
        display(plt)


        # y-column for currently engaged PMoA 
        c = Symbol("y_j_1_$(pmoa_idx)")

        # for every PMoA, verify that we have the effects where we need them
        for idx in 1:6
            # get the associated y-column
            c_idx = Symbol("y_j_1_$(pmoa_idx)")
            # if this is not the currently engaged pmoa
            if !(c_idx==c)
                # verify that y is close to 1 (== no effect)
                @test unique(0.99 .< sim[:,c_idx] .< 1.01) == [true]
            # if this is the currently engaged pmoa
            else
                # the highest concentration is equal to the median effective damage, 
                # so we expect y to bottom out around 0.5 for the given parameters (or max out at 2. for PMoA M)
                if pmoa_idx != 2
                    @test isapprox(0.5, minimum(sim[:,c]), atol = 0.01)
                else
                    @test isapprox(2.0, maximum(sim[:,c]), atol = 0.01)
                end
            end
        end

        # verify monotonically increasing effects on R - this should be the case for all PMoAs

        sim_end = @subset(sim, :t .== maximum(:t))
        sort!(sim_end, :R, rev = true)
        @test sort(sim_end.C_W_1) == sim_end.C_W_1
    end
end




#= 
code to check an individual PMoA by setting pmoa_idx.
this is not part of the tests. can be useful for debugging.

@testset "PMoA kappa" begin 
    global p = deepcopy(defaultparams)
<
    pmoa_idx = 6

    p.glb.t_max = 60.
    #p.glb.pathogen_inoculation_time = Inf

    p.glb.dX_in = [15., 15.]

    p.spc.H_j1 = 0.2
    p.spc.H_p = 40.

    p.spc.KD[pmoa_idx] = 1.
    p.spc.E[pmoa_idx] = 2. 
    p.spc.B[pmoa_idx] = 2.

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

 
end
=#