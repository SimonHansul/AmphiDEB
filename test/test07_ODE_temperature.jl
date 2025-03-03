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


using Plots 



AmphiDEB.defaultparams.spc.b_T


#@testset "Default parameters" begin
begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 80.
    p.glb.pathogen_inoculation_time = Inf
    p.glb.dX_in = 20.
    p.spc.H_p = 50.

    global sims = DataFrame()

    for T in [15., 20., 25.]
        p.glb.T = T + 273.15

        @time sim = AmphiDEB.ODE_simulator(
                p, 
                saveat = 1/240,
                alg = Tsit5()
                );

        sim[!,:E_mt_rel] = sim.E_mt ./ (sim.S + sim.E_mt)
        sim[!,:dI] = vcat(0, diff(sim.I))
        sim[!,:dI_rel] = sim.dI ./ (sim.S .^(2/3))
        sim[!,:W_tot] = sim.S .+ sim.E_mt 
        sim[!,:T] .= T

        append!(sims, sim)
    end
        
    @df sims plot(:t, :W_tot, group = :T, label = hcat(unique(:T)...), leg = true)
end














