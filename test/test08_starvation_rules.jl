using Pkg; Pkg.activate("test/")

using EcotoxSystems
using Plots, StatsPlots
using Test

using Revise
using AmphiDEB

p = deepcopy(AmphiDEB.defaultparams)


 
@testset "Reserve dynamics under starvation" begin
    p.glb.t_max = 365. * 1

    p.glb.dX_in = [1., 1.] # limit food availability
    p.spc.H_j1 = Inf # ignore metamorphosis
    p.spc.H_p = Inf # ignore adulthood
    
    @time sim = AmphiDEB.ODE_simulator(
        p
        )
    @df sim plot(
        plot(:t, [:S :E_mt], label = ["S" "E_mt"]),
        plot(:t, vcat(0, diff(:M)), label = ""),
        plot(:t, vcat(0, diff(:I)), label = ""), 
        xlabel = "t", ylabel = ["W" "dI" "dM"], 
        layout = (1,3), size = (1000,350)
    )

    # without additional starvation rules, 
    # change in reserve should go to 0. no endless reserve accumulation as individuals stop growing
    @test maximum(diff(sim.E_mt)[end-10:end]) <= 1e-4
            
end



