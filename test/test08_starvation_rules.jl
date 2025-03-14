using Pkg; Pkg.activate("test/")

using EcotoxSystems, AmphiDEB
using Plots, StatsPlots

p = deepcopy(AmphiDEB.defaultparams)

begin
    p.glb.t_max = 365. * 10
    p.glb.dX_in = [10., 10.]

    p.spc.H_j1 = Inf
    
    sim = AmphiDEB.ODE_simulator(p)
    @df sim plot(
        plot(:t, [:S :E_mt]),
        plot(:t, vcat(0, diff(:I))), 
        xlabel = "t", ylabel = ["W" "dI"]
    )
            
end

