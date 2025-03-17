
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
    gen_ind_params = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),
    kwargs...
    )

    EcotoxSystems.ODE_simulator(
        p;
        model = model,
        statevars_init = statevars_init,
        tstops = [p.glb.pathogen_inoculation_time, p.glb.medium_renewals...],
        gen_ind_params = gen_ind_params,
        callbacks = callbacks,
        kwargs...
    )
end

function IBM_simulator(
    # parameters
    p::ComponentVector;

    # default global model
    init_global_statevars = initialize_global_statevars,
    global_ode! = AmphiDEB_global!,
    global_rules! = AmphiDEB_global_rules!,
    
    # default individual model
    init_individual_statevars = initialize_individual_statevars,
    individual_ode! = AmphiDEB_individual!, 
    individual_rules! = default_individual_rules!,
    gen_ind_params = p -> EcotoxSystems.generate_individual_params(p; pth = p.pth),

    kwargs...
    )

    EcotoxSystems.IBM_simulator(
        p;

        # global model
        init_global_statevars = init_global_statevars,
        global_ode! = global_ode!, 
        global_rules! = global_rules!,
        
        # individual model
        individual_ode! = individual_ode!,
        init_individual_statevars = init_individual_statevars,
        individual_rules! = individual_rules!,
        gen_ind_params = gen_ind_params,

        kwargs...
    )

end