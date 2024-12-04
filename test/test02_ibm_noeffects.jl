using Pkg; Pkg.activate("AmphiDEB/test")

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
import EcotoxSystems: IBM_simulator


# Amphibian model without chemical or pathogen effects 
Amphibian_noeffects!(du, u, p, t) = begin

    DEBODE_global!(du, u, p, t)
    u.ind.y_j .= 1 # making sure that all effects are turned off
    u.ind.y_jP .= 1
    Amphibian_DEB!(du, u, p, t)

end

import AmphiDEB: IBM_simulator
import EcotoxSystems: default_individual_rules!

norm(x) = x ./ sum(x)

begin
    p = deepcopy(AmphiDEB.defaultparams)

    p.glb.t_max = 365*2
    p.glb.dX_in = 100.
    p.glb.k_V = 0.
    p.glb.N0 = 10

    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = 0 #365
    p.spc.h_S = 0.1

    @time sim = IBM_simulator(
        p; 
        individual_ode! = Amphibian_noeffects!, 
        init_individual_statevars = AmphiDEB.initialize_individual_statevars,
        showinfo = 30,  # update every 30 days
        saveat = 7, # saving weekly output
        dt = 1 # daily timestep - better to turn down to hourly for proper results
        )
    
    
    p_glb = @df sim.glb plot(:t, :N, xrotation = 45, xlabel = "Time [d]", ylabel = "N", title = "Abundance", leftmargin = 5mm)

    p_spc = @df sim.spc plot(
        groupedlineplot(:t, :S, :cohort, ylabel = "S", title = "Structural mass"), 
        groupedlineplot(:t, :cum_repro, :cohort, ylabel = "cR", title = "Cumulative reproduction"),
        groupedlineplot(:t, :f_X, :cohort, ylabel = "f(X)", ylim = (0.5, 1.01), title = "Scaled funct. response"), 
        groupedlineplot(:t, :H, :cohort, ylabel = "H", title = "Maturity"),
        xrotation = 45, xlabel = "Time [d]", titlefontsize = 10, 
        bottommargin = 5mm, topmargin = 5mm, leftmargin = 5mm, rightmargin = 5mm
    )
    
    plot(p_glb, p_spc, layout = grid(1,2, widths = norm([2/3, 1])), size = (1000,600))
end

