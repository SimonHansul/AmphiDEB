begin 
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
    import AmphiDEB: IBM_simulator

    import CSV
    norm(x) = x ./ sum(x)
end

begin
    p = deepcopy(defaultparams)

    p.glb.t_max = 365. * 2
    p.glb.dX_in = [100., 5000.]
    p.glb.k_V = [0., 0.]
    p.glb.N0 = 10

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = 0 #365
    p.spc.h_S = 0.1
    p.spc.H_p = 40.
    
    @time sim = AmphiDEB.IBM_simulator(
        p; 
        showinfo = 60,  # update every 30 days
        saveat = 7, # saving weekly output
        dt = 1 # daily timestep - better to turn down to hourly for proper results
    )

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
        #groupedlineplot(:t, [x[2] for x in :f_X], :cohort, ylabel = "f(X) terrestrial", title = "Scaled funct. response"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        groupedlineplot(:t, :A, :cohort, ylabel="A", title="Assimilation"),
        groupedlineplot(:t, :M, :cohort, ylabel="M", title="Maintenance cost"),
        maturity,
        life_stages,
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display
end

savefig("5jahre_dXin5000.png")