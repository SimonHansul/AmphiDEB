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

using AmphiDEB
norm(x) = x ./ sum(x)

function define_defaultparams()::EcotoxSystems.ComponentVector

    p = EcotoxSystems.ComponentVector(
    glb = AmphiDEB.defaultparams.glb, 
    pth = AmphiDEB.defaultparams.pth,
    spc = EcotoxSystems.ComponentVector(
        AmphiDEB.defaultparams.spc; 
        H_j1_prime = 1., 
        H_p_prime = 1.,
        watercontent_larvae = 0.93, 
        watercontent_juveniles = 0.85,
        time_since_birth = 15. # time since birth at the start of the experiment
    ))

    # setting global parameters

    p.glb.t_max = 365. 
    p.glb.pathogen_inoculation_time = Inf
    #p.glb.dX_in = 1e10 # ad libitum feeding conditions

    p.spc.Z = truncated(Normal(1, 0.17), 0, Inf)
    # propagation of zoom factor to H_j1 is turned off => we want variability in the transition to metamorphs
    p.spc.propagate_zoom.H_j1 = 0.

    p.spc.X_emb_int = 1. # ≈ initial dry mass of an egg (mg)


    p.spc.dI_max_lrv = 2.17
    p.spc.eta_AS_emb = 0.77
    p.spc.gamma = 0.85
    p.spc.H_j1 = 15.14
    p.spc.k_J_emb = 0.006
    p.spc.k_M_emb = 0.094
    p.spc.kappa_emb = 0.61


    p.spc.dI_max_juv = 3.28 
    p.spc.eta_AS_juv = 0.07
    p.spc.H_p = 2395.79
    p.spc.k_J_juv = 0.02
    p.spc.eta_AR = 0.16

    return p
end

@testset "Uninhibited growth" begin
    global p = define_defaultparams() 

    p.glb.t_max = 450.
    p.glb.dX_in = [1000., 1000.]
    p.glb.k_V = [0., 0.]
    p.glb.N0 = 10

    #p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = 0 #365
    p.spc.h_S = 0.1
    #p.spc.H_p = 40.
    
    @time sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60,  # update every 30 days
        saveat = 7, # saving weekly output
        dt = 1/24, # daily timestep - better to turn down to hourly for proper results
        record_individuals = false
        )
    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    
    CSV.write("sim_ind.csv", sim.spc)
    CSV.write("sim_glb.csv", sim.glb)

    #calculate the biomass of the population
    S_sum = zeros(Int(p.glb.t_max))
    timepoints = [i for i in 1:p.glb.t_max]

    j = 0
    for i in sim.spc.t
        j += 1
        S_sum[Int(floor(i)+1)] += sim.spc.S[j]
    end
    j = 0
    while j < size(S_sum)[1]
        j += 1
        if S_sum[j] == 0
            deleteat!(S_sum,j)
            deleteat!(timepoints, j)
            j -= 1
        end
    end

    #plot the food sources in one graph
    food = plot(sim.glb.t, sim.glb.X_1, xlabel="Time [d]", ylabel = "X", label="aquatic", title = "Food amount")
    plot!(sim.glb.t, sim.glb.X_2, label ="terrestrial")
    #plot!(sim.glb.t, [(i+1)*100 for i in sim.glb.t], label="dX_sum")
    plot!(legend=:topleft)

    #plot global state variables
    p_glb = @df sim.glb plot(
        plot(:t, :N, xrotation = 45, xlabel = "Time [d]", ylabel = "N", title = "Abundance"),
        #plot(timepoints, S_sum, xlabel = "Time [d]", ylabel = "total S", title = "biomass"),
        food,
        titlefontsize = 10, bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )

    #plot the life stages in one graph
    life_stages = groupedlineplot(sim.spc.t, sim.spc.embryo, sim.spc.cohort, xlabel="Time [d]", ylabel="population [%]", label="embryo", title="Life stages")
    groupedlineplot!(sim.spc.t,sim.spc.larva, sim.spc.cohort, label="larva")
    groupedlineplot!(sim.spc.t,sim.spc.metamorph, sim.spc.cohort, label="metamorph")
    groupedlineplot!(sim.spc.t,sim.spc.juvenile, sim.spc.cohort, label="juvenile")
    groupedlineplot!(sim.spc.t,sim.spc.adult, sim.spc.cohort, label="adult")
    #plot!(legend=:right)

    #plot maturity level with maturity thresholds in one graph
    maturity = groupedlineplot(sim.spc.t, sim.spc.H, sim.spc.cohort, xlabel="Time [d]", ylabel = "H", title = "Maturity")
    plot!(sim.spc.t, [p.spc.H_j1*i/i for i in sim.spc.t])
    plot!(sim.spc.t, [p.spc.H_p*i/i for i in sim.spc.t])

    #plot individual state variables
    p_spc = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        #groupedlineplot(:t, [x[1] for x in :f_X], :cohort, ylabel = "f(X) aquatic", title = "Scaled funct. response"), 
        groupedlineplot(:t, [x[2] for x in :f_X], :cohort, ylabel = "f(X) terrestrial", title = "Scaled funct. response"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        groupedlineplot(:t, :A, :cohort, ylabel="A", title="Assimilation"),
        groupedlineplot(:t, :M, :cohort, ylabel="M", title="Maintenance cost"),
        maturity,
        #life_stages,
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display

    # check expected order of magnitude of population size after 1.5 years of simulation with defaultparams
    @test 500 < maximum(sim.glb.N) < 1000
end

savefig("discoglossus_test.png")
#p.glb.t_max = 300
#VSCodeServer.@profview_allocs AmphiDEB.IBM_simulator(
#    p; 
#    showinfo = 60,  # update every 30 days
#    saveat = 7, # saving weekly output
#    dt = 1/24 # daily timestep - better to turn down to hourly for proper results
#    )

@testset "Simulation with density-dependence" begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365. * 3
    p.glb.dX_in = [500., 500.]
    p.glb.k_V = [0.1, 0.1]
    p.glb.N0 = 100

    p.spc.X_emb_int = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(365., 36.5), 0, Inf)
    p.spc.H_p = 50.
    p.spc.K_X_lrv = p.spc.K_X_juv = 20.
    p.spc.h_b = -log(1 - 0.005)

    p.spc.h_S = -log(0.75)
    p.spc.S_rel_crit = 0.5

    @time global sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60, 
        saveat = 1, 
        dt = 1/24, 
        record_individuals = false 
        )

    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years

    p1 = @df sim.glb plot(
        :t, :N, 
        xrotation = 45, 
        xlabel = "Time [d]", ylabel = "N", title = "Total abundance", 
        leftmargin = 5mm, bottommargin = 5mm
        )

    p2 = @df sim.glb plot(
        :t, [:N_emb :N_lrv :N_mt :N_juv :N_ad], 
        label = ["Embryos" "Larvae" "Metamorphs" "Juveniles" "Adults"], 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance per life stage", 
        leg = :outertopright, xrotation = 45, bottommargin = 5mm
    )
    
    plot(p1, p2, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display
end
