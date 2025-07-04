# modelcollection.jl
# provides an overview of pre-defined models as dict

ODE_models = Dict(
    :complete => Dict(
        :M1_complete_with_linear_TD => M1_complete_ODE_with_linear_TD!, 
        :M1_complete_ODE_with_loglogistic_TD => M1_complete_ODE_with_loglogistic_TD!, 
        :M2_complete_ODE_with_loglogistic_TD => M2_complete_ODE_with_loglogistic_TD!
    )
)