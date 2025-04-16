using Pkg; Pkg.activate("test")
using Test
using Distributions
using OrdinaryDiffEq

using Plots, StatsPlots, Plots.Measures
default(leg = false, lw = 1.5)

include("testutils.jl")

using DataFrames, DataFramesMeta, CSV
using StatsBase

using EcotoxSystems
import EcotoxSystems: constrmvec

using Revise

using AmphiDEB
@time import AmphiDEB: defaultparams, ODE_simulator, Amphibian_DEB!, AmphiDEB_ODE!


include("scripts/timevar_test_utils.jl")


exposure = CSV.read("test/data/timevar_test/exposure.csv", DataFrame, types = Float64)



@df exposure plot(
    plot(:t, :C_W_1, group = :scenario, leg = true, label = hcat(unique(:scenario)...), legendtitle = "Scenario", legendtitlefontsize = 8),
    plot(:t, :C_W_2, group = :scenario, leg = true, label = hcat(unique(:scenario)...), legendtitle = "Scenario", legendtitlefontsize = 8)
)

df = @subset(exposure, :scenario .== 1)


interp = linear_interpolation(df.t, df.C_W_1)

interp(14)


p = deepcopy(AmphiDEB.defaultparams)

exposure_df = DataFrame(
    scenario = [1, 1, 1, 1],
    t = [0, 5, 5, 10],
    C_W = [0, 0, 30., 30.]
)

p.glb.C_W = [exposure_df]


exposure_df = DataFrame()


p.glb.C_W