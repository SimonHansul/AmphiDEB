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
norm(x) = x ./ sum(x)

# 
# TODO: implement multiple food resources 
#   - Adjust global state variables
#       - X is a vector
#       - optinal: C_W is vector
#   - Adjust individual state variables
#       - f_X is a vector
#       - I_p is a vector
#       - in theory, we also have to change shape of D_z etc...ignoring for now
#       
#   - Adjust global parameters
#       - dX_in is a vector
#       - k_V is a vector 
#       - V_patch is a vector OR add A_patch
#   - Adjust species parameters:
#       - K_X_lrv, K_Xjuv, dI_max_lrv etc. can remain as they are
#       - assume that eta_IA remains constant
#   
#   - Adjust derivatives
#       - Line 89: f_X needs to use correct X, K_X, and V_patch 
#       - Lines 91-94: Make sure that corect f_X value is applied, depending on the life stage
#   - Check that zoom factor is still applied correctly

# TODO: later (optional)
#   - add seasonal temperature change (e.g. sine curve)
#   - difference between water and air temperature?


#warum ist v_patch = [0.0,0.0]? ist weg?
#ecotox.ibmschedule model_step X = max(0,X) element wise machen 

#schlechtere maintenance und assimilation verhindern erreichen des maturity threshold

#todo 6.2.:
#dose repsonse plots nochmal mit niedrigerem maturity threshold machen
#tests zur 2d nachrung machen


##25.02.
#mit Simon treffen. woher Werte? realistische bereiche? aktuelle ergebnisse
#pathogen und tktd fixen
#dosis wirkungskurven für einzelne generationen
#nahrungsverfügbarkeit deutlich erhöhen
#hintergrundmortalität