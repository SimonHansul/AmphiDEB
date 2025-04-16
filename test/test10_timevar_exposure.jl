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

exposure_scenarios = CSV.read("test/data/timevar_test/exposure.csv", DataFrame, types = Float64)

@df exposure_scenarios plot(
    plot(:t, :C_W_1, group = :scenario, xlim = (0,14), leg = true, label = hcat(unique(:scenario)...), legendtitle = "Scenario", legendtitlefontsize = 8),
    plot(:t, :C_W_2, group = :scenario, xlim = (0,14), leg = true, label = hcat(unique(:scenario)...), legendtitle = "Scenario", legendtitlefontsize = 8), 
    xlabel = "Time (d)", ylabel = ["C_W" ""], title = ["Stressor 1" "Stressor 2"]
)


using Interpolations

#function generate_interpolators(
#    exposure_scenarios::AbstractDataFrame,
#    t_max::Real, 
#    interpolation_kind::AbstractString = "linear",
#    plot_interpolations::Bool = true
#    )

t_max = 14.

exposure_cols = filter(x -> occursin("C_W_", x), names(exposure_scenarios))

i = 1
scenario = 4

df = @subset(exposure_scenarios, :scenario .== scenario)
sort!(exposure_scenarios, :t)
df.t = Interpolations.deduplicate_knots!(df.t)

interp_1 = linear_interpolation(df.t, df[:,exposure_cols[1]])
interp_2 = linear_interpolation(df.t, df[:,exposure_cols[2]])

p.glb.C_W[1] = interp_1
p.glb.C_W[2] = interp_2

plot(
    plot(x -> interp_1(x), xlim = (0,t_max), xlabel = "t", ylabel = exposure_cols[i], title = "Scenario: $(scenario)"),
    plot(x -> interp_2(x), xlim = (0,t_max), xlabel = "t", ylabel = exposure_cols[i], title = "Scenario: $(scenario)")
)
@df df scatter!(:t, :C_W_1)
@df df scatter!(:t, :C_W_2)


sim = AmphiDEB.ODE_simulator(p)

EcotoxSystems.generate_individual_params(p)

ind = EcotoxSystems.getval.(p.spc) |> EcotoxSystems.propagate_zoom


p.glb.C_W

EcotoxSystems.ComponentVector(
    glb = p.glb
)


@test "Running a single time-variable exposure scenario" begin


end

interp(10.)


p = deepcopy(AmphiDEB.defaultparams)

exposure_df = DataFrame(
    scenario = [1, 1, 1, 1],
    t = [0, 5, 5, 10],
    C_W = [0, 0, 30., 30.]
)

p.glb.C_W = [exposure_df]


exposure_df = DataFrame()


p.glb.C_W