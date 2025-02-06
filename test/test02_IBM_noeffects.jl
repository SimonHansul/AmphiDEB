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

using Revise

@time import AmphiDEB: defaultparams, ODE_simulator, Amphibian_DEB!, AmphiDEB_ODE!
using AmphiDEB
norm(x) = x ./ sum(x)

@testset "Uninhibited growth" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 450.
    p.glb.dX_in = 1000.
    p.glb.k_V = 0.
    p.glb.N0 = 10

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = 1
    p.spc.h_S = 0.
    p.spc.H_p = 50.

    @time global sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60,  # update every 30 days
        saveat = 7, # saving weekly output
        dt = 1/24 # daily timestep - better to turn down to hourly for proper results
        )
    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    p_glb = @df sim.glb plot(
        :t, :N, 
        xrotation = 45, 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance", 
        leftmargin = 5mm
        )

    p_spc = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X)", ylim = (0.5, 1.01), title = "Scaled funct. response"), 
        groupedlineplot(:t, :H, :cohort, ylabel = "H", title = "Maturity"),
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display


    @test 2e4 < maximum(sim.glb) < 2e5
end


# FIXME: when we have a yearly breeding pulse, the second generation just stays small and keeps hanging out...

@testset "Simulation with density-dependence" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 10. #365. * 3
    p.glb.dX_in = 500.
    p.glb.k_V = 0.1
    p.glb.N0 = 10

    p.spc.X_emb_int = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = 365.
    p.spc.H_p = 50.
    p.spc.K_X_lrv = p.spc.K_X_juv = 20.

    p.spc.h_S = -log(0.5)
    p.spc.S_rel_crit = 0.5

    @time global sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60, # print update every so many days 
        saveat = 1, # saving weekly output
        dt = 1/24 # daily timestep - better to turn down to hourly for proper results
        )
    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    p_glb = @df sim.glb plot(
        :t, :N, 
        xrotation = 45, 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance", 
        leftmargin = 5mm
        )

    p_spc = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        groupedlineplot(:t, :X_emb, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X)", title = "Scaled funct. response"), 
        groupedlineplot(:t, :H, :cohort, ylabel = "H", title = "Maturity"),
        rotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display
end


plot(
    (@df sim.glb plot(:t, :X ./ p.glb.V_patch)),
    (@df sim.spc plot(:t, :f_X))
)

sim.spc[!,:lifestage] = [argmax([r.embryo, r.larva, r.metamorph, r.juvenile, r.adult]) for r in eachrow(sim.spc)]

popstruct = combine(groupby(sim.spc, :t)) do df
    DataFrame(
        num_emb = sum(df.lifestage .== 1), 
        num_larva = sum(df.lifestage .== 2),
        num_metamorph = sum(df.lifestage .== 3),
        num_juvenile = sum(df.lifestage .== 4),
        num_adult = sum(df.lifestage .== 5),
    )
end

@df popstruct plot(
    :t, [:num_emb :num_larva :num_metamorph :num_juvenile :num_adult], 
    label = ["embryos" "larvae" "metamorphs" "juveniles" "adults"]
    )

