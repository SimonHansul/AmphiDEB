using Pkg; Pkg.activate("profiling")

using BenchmarkTools
using Revise
using AmphiDEB

#=
## log
- Initial benchmark 
    - median = 12.5 ms (9-48)
    - 17.9 MB allocs 
- added @inbounds to growth derivative
    - 11.4 ms
    - 14.6 MB
- added pre-indexing in determine_life_stage
    - 9.5 ms
    - 12.3 MB
- added inbounds to global deriv
    - 9.3 ms
    - 11.2 MB
- added pre-indexing to life stage effects
    - 9.2 ms
    - 10 MB
- added @views to TKTD
    - 10.7 MB
    - --> removed again
- changed return type of LL2 from Real to Float64
    - no change
- same for other LL functions
    - 10.5 MB
- removing dispatch from LL2(x, p1, p2) to LL2(x, p::Tuple)
    - no change
- removing call to complex from LL2
     - 8.4 ms
     - 10.7 MB
- removing call to complex from LL2 + removing ifelse
    - no considerable changes...
- applying inbounds to entire M1_ingestion function
    - no change
- added usual optims to pathogen function
     - no improvement
=#

begin
    p = deepcopy(AmphiDEB.defaultparams)
    AmphiDEB.ODE_simulator(p);
    @benchmark AmphiDEB.ODE_simulator(p)
end

#=
IBM does barely more allocs than the ODE

in global rules,
replaced map with generateor 
- massive imporvement!
=#

begin
    AmphiDEB.IBM_simulator(p)
 
    @benchmark AmphiDEB.IBM_simulator(p)
end

# simulation with map() --> 417s, 170 GB 
# simulation with generator --> same
# simulation without either --> minimal change
#

begin
    global p = deepcopy(defaultparams)

    p.glb.t_max = 365. * 3
    p.glb.dX_in = [100., 5000.]
    p.glb.k_V = [0.1, 0.1]
    p.glb.N0 = 100

    p.spc.X_emb_int = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.Z = truncated(Normal(1, 0.1), 0, Inf)
    p.spc.tau_R = truncated(Normal(365., 36.5), 0, Inf)
    p.spc.H_p = 50.
    p.spc.K_X_lrv = p.spc.K_X_juv = 20.
    p.spc.h_b = -log(1 - 0.005)

    p.spc.h_S = -log(0.75)
    p.spc.S_rel_crit = 0.5

    @time sim = AmphiDEB.IBM_simulator(p; showinfo = 60, record_individuals = false);

    nothing
end

@df sim.glb plot(:t, :N)

@profview_allocs AmphiDEB.IBM_simulator(p; record_individuals = false)
