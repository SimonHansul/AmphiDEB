# derivatives.jl
# model functions which are used by all models

function Pathogen_growth!(du, u, p, t)::Nothing
    
    du.glb.P_Z = - p.pth.mu * u.glb.P_Z # change in zoospores in environment
    
    return nothing
end


function Pathogen_Infection!(du, u, p, t)::Nothing

    # life-stage specificity of the zoospore-host encounter rate 
    # - currently turned off to simplify things 
    # this needs to be re-evaluated. tadpoles do get infected but are unlikely to die 
    # maybe only make the effect parameters life stage-specific?
    gamma = p.pth.gamma 

    # growth and killing of sporangia 
    du.ind.P_S = p.pth.v0 * gamma * u.glb.P_Z + p.pth.v0 * p.pth.eta * p.pth.f * u.ind.P_S - (p.pth.sigma0 + p.pth.sigma1 * p.ind.Chi * u.ind.P_S) * u.ind.P_S
    
    # feedback with zoospore population 
    du.glb.P_Z += p.pth.eta * (1-p.pth.f) * u.ind.P_S # release of spores
    du.glb.P_Z -= gamma * u.glb.P_Z # encystment

    # relative response to pathogen 
    try
    
        @. u.ind.y_jP = EcotoxSystems.LL2(u.ind.P_S/(Complex(u.ind.S^(2/3)).re), p.ind.e_P, p.ind.b_P)

    catch
        println((u.ind.S))
        error()
    end

    u.ind.y_jP[2] /= u.ind.y_jP[2]^2 # converting a monotonically decreasing to increasing response

    return nothing
end

# starting with the definition of callbacks to trigger life stage transitions and global events (e.g. pathogen inoculation)
# definition of the callback to incoulate pathogen at a given time-point

condition_inoculation(u, t, integrator) = integrator.p.glb.pathogen_inoculation_time - t
function effect_inoculation!(integrator) 
    integrator.u.glb.P_Z = integrator.p.glb.pathogen_inoculation_dose
end

# definition of the callback for media renewals
condition_renewal(u, t, integrator) = prod(integrator.p.glb.medium_renewals  .- t)  
function effect_renewal!(integrator)
    integrator.u.glb.P_Z = 0. # when renewal occurs, set zoospores to 0
end


"""
    ODEcallbacks()

Callbacks to use in the ODE system
"""
function AmphODE_callbacks()

    # life-stage callbacks are discontinued for this model atm, 
    # I added the definition of life stages directly to the derivatives
    # should also consider to define the life stage tranistions as smooth sigmoid switches, at least for metamorphs

    #cb_larva = ContinuousCallback(condition_larva, effect_larva!)
    #cb_metamorph = ContinuousCallback(condition_metamorph, effect_metamorph!)
    #cb_juvenile = ContinuousCallback(condition_juvenile, effect_juvenile!)
    #cb_adult = ContinuousCallback(condition_adult, effect_adult!)

    cb_inoculation = ContinuousCallback(condition_inoculation, effect_inoculation!)
    cb_renewal = ContinuousCallback(condition_renewal, effect_renewal!)

    return CallbackSet(
        cb_inoculation, 
        cb_renewal
        )

end

function determine_life_stage!(du, u, p, t)::Nothing

    #TODO: between larvae and metamorphs, there should be a smooth transition using a sgmoid function depending on H
    
    u.embryo = sig(u.X_emb, 0, 0, 1, beta = 1e3) #u.X_emb > 0 # if there is still vitellus left, we are in embryo stage
    u.larva = sig(u.X_emb, 0, 1, 0, beta = 1e3) * sig(u.H, p.H_j1 * u.y_j[5], 1, 0) #(u.H <= p.H_j1 * u.y_j[5]) # if the embryo is used up but the next maturity threshold is not reached, we are in larval stage
    u.metamorph =  sig(u.H, p.H_j1 * u.y_j[5], 0, 1) * sig(u.E_mt, 0, 0, 1, beta = 1e3) # above the maturity threshold for Gosner stage 42, while there is still metamorphic reserve left, we are in metamorph stage
    u.juvenile = sig(u.H, p.H_j1 * u.y_j[5], 0, 1) * sig(u.E_mt, 0, 1, 0, beta = 1e3) * sig(u.H, p.H_p, 1, 0) # after metamorphosis but below the threshold for puberty, we are in juvenile stage
    u.adult = sig(u.H, p.H_p, 0, 1)

    return nothing
end


# temperature correction
function y_T!(du, u, p, t)::Nothing
    
    u.ind.y_T = exp((p.ind.T_A / p.ind.T_ref) - (p.ind.T_A / p.glb.T)) 

    return nothing
end


"""
    ingestion!(du, u, p, t)::Nothing

Life stage-specific calculation of ingestion rate for amphibians. <br>
Note that this function requires all components to be included in the arguments `du`, `u` and `p`, so that we can access the external food concentration. <br>
Hence, we have `u.ind`, `u.glb`, `p.ind`, etc., instead of simply `u` and `p`.
"""
function ingestion!(du, u, p, t)::Nothing

    K_X = ((u.ind.larva + u.ind.metamorph) * p.ind.K_X_lrv) + ((u.ind.juvenile + u.ind.adult) * p.ind.K_X_juv)
    u.ind.f_X = (u.glb.X / p.glb.V_patch) / ((u.glb.X / p.glb.V_patch) + K_X)

    dI_emb = u.ind.embryo * (Complex(u.ind.S)^(2/3)).re * p.ind.dI_max_emb * u.ind.y_T
    dI_mt = u.ind.metamorph * u.ind.f_X * p.ind.dI_max_lrv * (u.ind.E_mt / u.ind.E_mt_max) * (Complex(u.ind.S)^(2/3)).re * u.ind.y_T
    dI_lrv = u.ind.larva * u.ind.f_X * p.ind.dI_max_lrv * (Complex(u.ind.S)^(2/3)).re * u.ind.y_T
    dI_juv_ad = (u.ind.juvenile + u.ind.adult) * u.ind.f_X * p.ind.dI_max_juv * (Complex(u.ind.S)^(2/3)).re * u.ind.y_T
    
    du.ind.I = dI_emb + dI_mt + dI_lrv + dI_juv_ad
    du.ind.X_emb = -dI_emb
    
    du.glb.X -= du.ind.I

    # assimilation flux
    du.ind.A = du.ind.I * p.ind.eta_IA * u.ind.y_j[3] * u.ind.y_jP[3]

    return nothing 
end

"""
    maintenance!(du, u, p, t)::Nothing 

Calculation of maintenance fluxes for amphibians.
"""
function maintenance!(du, u, p, t)::Nothing 

    k_M = (u.embryo + u.larva + u.metamorph) * p.k_M_emb + (u.juvenile + u.adult) * p.k_M_juv
    k_J = (u.embryo + u.larva + u.metamorph) * p.k_J_emb + (u.juvenile + u.adult) * p.k_J_juv

    du.M = u.S * k_M * u.y_j[2] * u.y_jP[2] * u.y_T # somatic maintenance flux
    du.J = u.H * k_J * u.y_j[2] * u.y_jP[2] * u.y_T # maturity maintenance flux

    return nothing
end

function reproduction!(du, u, p, t, kappa)::Nothing
    
    du.R = u.adult * clipneg(p.eta_AR * u.y_j[4] * u.y_jP[4] * ((1 - kappa) * du.A - du.J))  # reproduction flux

    return nothing
end
