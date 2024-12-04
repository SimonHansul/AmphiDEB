using Pkg; Pkg.activate("AmphiDEB/test")

using Test
using Plots, StatsPlots, Plots.Measures
using Distributions
using OrdinaryDiffEq
default(leg = false, lw = 1.5)

using Chain
using DataFrames, DataFramesMeta
using StatsBase

using Revise


include("test01_base_model_noeffects.jl")
include("test02_pathogens.jl")