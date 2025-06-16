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
        :t, [:N :N_emb :N_lrv :N_mt :N_juv :N_ad], 
        label = ["Abundance" "Embryos" "Larvae" "Metamorphs" "Juveniles" "Adults"], 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance per life stage", 
        leg = :topright, xrotation = 45, bottommargin = 5mm
    )
    #groupedlineplot(sim.spc.t, sim.spc.embryo, sim.spc.cohort, xlabel="Time [d]", ylabel="population [%]", label="embryo", title="Life stages")
    #groupedlineplot!(sim.spc.t,sim.spc.larva, sim.spc.cohort, label="larva")
    #groupedlineplot!(sim.spc.t,sim.spc.metamorph, sim.spc.cohort, label="metamorph")
    #groupedlineplot!(sim.spc.t,sim.spc.juvenile, sim.spc.cohort, label="juvenile")
    #groupedlineplot!(sim.spc.t,sim.spc.adult, sim.spc.cohort, label="adult")
    #plot!(legend=:right)

    #plot maturity level with maturity thresholds in one graph
    maturity = groupedlineplot(sim.spc.t, sim.spc.H, sim.spc.cohort, xlabel="Time [d]", ylabel = "H", title = "Maturity")
    plot!(sim.spc.t, [p.spc.H_j1*i/i for i in sim.spc.t])
    plot!(sim.spc.t, [p.spc.H_p*i/i for i in sim.spc.t])

    #plot individual state variables
    p_spc = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        #groupedlineplot(:t, [x[1] for x in :f_X], :cohort, ylabel = "f(X) aquatic", title = "Scaled funct. response"), 
        #groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X) terrestrial", title = "Scaled funct. response"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        #groupedlineplot(:t, :A, :cohort, ylabel="A", title="Assimilation"),
        #groupedlineplot(:t, :R, :cohort, ylabel="R", title="Reproduction buffer"),
        maturity,
        food,
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display

    # check expected order of magnitude of population size after 1.5 years of simulation with defaultparams
    @test 100 < maximum(sim.glb.N) < 1000
end

savefig("discoglossus_test.png")
#@testset "Simulation with density-dependence" 
begin
    global p = define_defaultparams() #deepcopy(defaultparams)

    p.glb.t_max = 365. * 10
    p.glb.dX_in = [5000., 5000.]
    p.glb.k_V = [0.1, 0.1]
    p.glb.N0 = 100

    #p.spc.X_emb_int = truncated(Normal(1., 0.1), 0, Inf)
    #p.spc.Z = truncated(Normal(1., 0.1), 0, Inf)
    #p.spc.tau_R = 0.0#truncated(Normal(365., 36.5), 0, Inf)
    #p.spc.H_p = 50.
    #p.spc.K_X_lrv = p.spc.K_X_juv = 20.
    p.spc.h_b = -log(1 - 0.001)

    p.spc.h_S = -log(0.75)
    p.spc.S_rel_crit = 0.5

    #modes of action
    p.spc.eta_IA *= 1.025 
    #p.spc.k_M_emb *= 1.05
    #p.spc.k_M_juv *= 1.05
    #p.spc.k_J_emb *= 1.05
    #p.spc.k_J_juv *= 1.05
    #p.spc.H_j1 *= 0.9

    @time global sim = @replicates AmphiDEB.IBM_simulator(
        p; 
        showinfo = 365, 
        saveat = 1, 
        dt = 1/24, 
        record_individuals = true 
    ) 10
    
    #CSV.write("sim_ind.csv", sim.spc)
    #CSV.write("sim_glb.csv", sim.glb)
    repetitions = maximum(sim.glb.replicate)

    maxi = []
    mini = []
    avg = []
    sigma = []

    mean_embryotime = zeros(repetitions)
    mean_larvaltime = zeros(repetitions)
    mean_metamorphtime = zeros(repetitions)
    mean_juveniletime = zeros(repetitions)
    mean_adulttime = zeros(repetitions)

    max_embryotime = zeros(repetitions)
    max_larvaltime = zeros(repetitions)
    max_metamorphtime = zeros(repetitions)
    max_juveniletime = zeros(repetitions)
    max_adulttime = zeros(repetitions)
    
    replic_glb = groupby(sim.glb, :replicate)
    population = combine(replic_glb, :N => mean => :mean, :N => std => :std, :N => maximum => :max, :N => minimum => :min) 
    avg = @. round(population[!,:mean],digits=2)
    sigma = @. round(population[!,:std],digits=2)
    maxi = @. round(population[!,:max],digits=2)
    mini = @. round(population[!,:min],digits=2)
    #abundance
    println("Maximum: $(maxi)")
    println("Maximum: $(round(mean(maxi),digits=2))")
    println("Minimum: $(mini)")
    println("Minimum: $(round(mean(mini),digits=2))")
    println("Mean: $(avg)")
    println("Mean: $(round(mean(avg),digits=2))")
    println("Std: $(sigma)")
    println("Std: $(round(mean(sigma),digits=2))")
    
    #@test 2500 <= sim.glb.N[end] <= 3500 # expected abundance after 2 years
    
    #N_replicates = @df sim.glb lineplot(:t, :N, xlabel="time [d]", ylabel="abundance [-]", label="m_thr -10%")
    #=
    p2 = @df sim.glb plot(
        :t, [:N :N_emb :N_lrv :N_mt :N_juv :N_ad],
        label = ["Abundance" "Embryos" "Larvae" "Metamorphs" "Juveniles" "Adults"], 
        xlabel = "Time [d]", ylabel = "N", title = "Abundance per life stage", 
        leg = :topright, xrotation = 45, bottommargin = 5mm
    )
    
    #plot maturity level with maturity thresholds in one graph
    maturity = groupedlineplot(sim.spc.t, sim.spc.H, sim.spc.cohort, xlabel="Time [d]", ylabel = "H", title = "Maturity")
    plot!(sim.spc.t, [p.spc.H_j1*i/i for i in sim.spc.t])
    plot!(sim.spc.t, [p.spc.H_p*i/i for i in sim.spc.t])

    #plot the food sources in one graph
    food = plot(sim.glb.t, sim.glb.X_1, xlabel="Time [d]", ylabel = "X", label="aquatic", title = "Food amount")
    plot!(sim.glb.t, sim.glb.X_2, label ="terrestrial")
    plot!(legend=:right)
    
    #plot individual state variables
    p3 = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        #groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X)", title = "Scaled funct. response"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        #groupedlineplot(:t, :A, :cohort, ylabel="A", title="Assimilation"),
        #groupedlineplot(:t, :R, :cohort, ylabel="R", title="Reproduction buffer"),
        maturity,
        food,
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p2, p3, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600)) |> display
    =#

    #life stages
    replic_spc = groupby(sim.spc, :replicate)
    for i in 1:repetitions
        individuals = groupby(replic_spc[(replicate=i,)], :id)
        life_stages = combine(individuals, [:embryo, :larva, :metamorph, :juvenile, :adult] => ((e,l,m,j,a) -> (a=sum(e), b=sum(l), c=sum(m), d=sum(j), e=sum(a))) => AsTable)
        mean_embryotime[i] = round(mean(filter(x-> x!=0, life_stages[!,"a"])),digits=2)
        max_embryotime[i] = round(maximum(life_stages[!,"a"]), digits=2)
        mean_larvaltime[i] = round(mean(filter(x-> x!=0, life_stages[!,"b"])),digits=2)
        max_larvaltime[i] = round(maximum(life_stages[!,"b"]), digits=2)
        mean_metamorphtime[i] = round(mean(filter(x-> x!=0, life_stages[!,"c"])),digits=2)
        max_metamorphtime[i] = round(maximum(life_stages[!,"c"]), digits=2)
        mean_juveniletime[i] = round(mean(filter(x-> x!=0, life_stages[!,"d"])),digits=2)
        max_juveniletime[i] = round(maximum(life_stages[!,"d"]), digits=2)
        mean_adulttime[i] = round(mean(filter(x-> x!=0, life_stages[!,"e"])),digits=2)
        max_adulttime[i] = round(maximum(life_stages[!,"e"]), digits=2)
    end
    println("Embryo_mean t=$(mean_embryotime)")
    println("Embryo_max t=$(max_embryotime)")

    println("larval_mean t=$(mean_larvaltime)")
    println("larval_max t=$(max_larvaltime)")

    println("metamorph_mean t=$(mean_metamorphtime)")
    println("metamorph_max t=$(max_metamorphtime)")

    println("juvenile_mean t=$(mean_juveniletime)")
    println("juvenile_max t=$(max_juveniletime)")

    println("adult_mean t=$(mean_adulttime)")
    println("adult_max t=$(max_adulttime)")

    #plot!()
end
plot!()
#savefig("discoglossus_100replicates.png")
#N_replicates = @df sim.glb lineplot(:t, :N, xlabel="time [d]", ylabel="abundance [-]", label="control", color="grey")
plot(N_replicates)
lineplot!(sim.glb.t, sim.glb.N, label="assimilation +2.5%")
plot!(legend=:topleft)
#savefig("discoglossus_control_assimilation2.png")

S_max_lrv = ((p.spc.kappa_emb*p.spc.dI_max_lrv*p.spc.eta_IA)/p.spc.k_M_emb)^3
S_max_ad = ((p.spc.kappa_juv*p.spc.dI_max_juv*p.spc.eta_IA)/p.spc.k_M_juv)^3
H_eq_lrv = ((1-p.spc.kappa_emb)*p.spc.dI_max_lrv*p.spc.eta_IA*S_max_lrv^(2/3))/p.spc.k_J_emb #4404.309199105734
H_eq_ad = ((1-p.spc.kappa_juv)*p.spc.dI_max_juv*p.spc.eta_IA*S_max_ad^(2/3))/p.spc.k_J_juv #2938.9853028115035
p.spc.H_j1/H_eq_lrv #metamorph_thr = 0.003437542487497036
p.spc.H_p/H_eq_ad #puberty_thr = 0.8151759036386231
