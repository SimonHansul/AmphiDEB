
"""
A wrapper around the `EcotoxSystems.jl` ODE_simulator 
for application with the AmphiDEB model.
"""
function ODE_simulator(
    p::ComponentVector; 
    model = AmphiDEB_ODE_M1!, 
    callbacks = AmphODE_callbacks(), 
    statevars_init = initialize_statevars,
    kwargs...
    )

    EcotoxSystems.ODE_simulator(
        p;
        model = model,
        statevars_init = statevars_init,
        tstops = [p.glb.pathogen_inoculation_time, p.glb.medium_renewals...],
        ind_params_init = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),
        callbacks = callbacks, 
        kwargs...
    )
end


function IBM_simulator(
    p::ComponentVector;
    individual_ode! = AmphiDEB_ODE_M1!, 
    init_individual_statevars = initalize_individual_statevars,
    callbacks = AmphODE_callbacks(), 
    statevars_init = initialize_statevars,
    kwargs...
    )

end