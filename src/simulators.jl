
"""
    ODE_simulator(
        p::ComponentVector; 
        model = AmphiDEB_ODE_M1!, 
        callbacks = AmphODE_callbacks(), 
        statevars_init = initialize_statevars,
        kwargs...
        )
A wrapper around the `EcotoxSystems.jl` ODE_simulator 
for application with the AmphiDEB model.
"""
function ODE_simulator(
    p::ComponentVector; 
    model = AmphiDEB_ODE!, 
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
    individual_ode! = AmphiDEB_ODE!, 
    individual_rules! = default_individual_rules!,
    init_global_statevars = initialize_global_statevars,
    init_individual_statevars = initialize_individual_statevars,
    global_rules! = EcotoxSystems.default_global_rules!,
    global_ode! = EcotoxSystems.DEBODE_global!,
    kwargs...
    )

    EcotoxSystems.IBM_simulator(
        p;
        individual_ode! = individual_ode!,
        init_global_statevars = init_global_statevars,
        init_individual_statevars = init_individual_statevars,
        individual_rules! = individual_rules!,
        global_ode! = global_ode!, 
        global_rules! = global_rules!
    )

end