

const X_EMB_INT_REL = 1e-3


function initialize_individual_statevars(p::ComponentVector; kwargs...)::ComponentVector
    ComponentVector(
        EcotoxSystems.initialize_individual_statevars(p); # state variables of the default model
        larva = 0, # additional life stages
        metamorph = 0,
        E_mt = 1e-10, # metamorphic reserve
        E_mt_max = 1e-10, # maximum metamorphic reserve
        P_S = 0, # pathogen sporangia 
        # this should use EcotoxSystems.constrmvec to construct a mutable static array
        # it doesn't work right now, not sure why, so I use regular arrays 
        # maybe slows down simulation a little compared to static array, but otherwise the same
        y_jP = ones(4), #EcotoxSystems.constrmvec(Vector(p.ind.e_P), fillval = 1) # response to pathogen per PMoA
        kwargs...
    )
end

initialize_pathogen_statevars(p) = ComponentVector(P_Z = 0)

function initialize_statevars(p::ComponentVector)

    global_statevars = ComponentVector(EcotoxSystems.initialize_global_statevars(p))
    pathogen_statevars = initialize_pathogen_statevars(p)
    individual_statevars = initialize_individual_statevars(p)
    
    return ComponentVector(
        glb = global_statevars,
        pth = pathogen_statevars,
        ind = individual_statevars
    )
end