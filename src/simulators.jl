
"""
A wrapper around the `EcotoxSystems.jl` ODE_simulator.
"""
function ODE_simulator(
    p::ComponentVector; 
    model = AmphiDEB_ODE_M1!, 
    callbacks = AmphODE_callbacks(), 
    kwargs...
    )

    EcotoxSystems.ODE_simulator(
        p;
        model = model,
        statevars_init = initialize_statevars,
        tstops = [p.glb.pathogen_inoculation_time, p.glb.medium_renewals...],
        ind_params_init = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),
        callbacks = callbacks, 
        kwargs...
    )
end
