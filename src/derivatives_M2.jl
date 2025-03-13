# derivatives_M2.jl
# derivatives for an alternative model version where E_mt is viewed as a sub-compartment of structure
"""
    ingestion!(du, u, p, t)::Nothing

Life stage-specific calculation of ingestion rate for amphibians. <br>
Note that this function requires all components to be included in the arguments `du`, `u` and `p`, so that we can access the external food concentration. <br>
Hence, we have `u.ind`, `u.glb`, `p.ind`, etc., instead of simply `u` and `p`.
"""
function ingestion_M2!(du, u, p, t)::Nothing

    K_X = ((u.ind.larva + u.ind.metamorph) * p.ind.K_X_lrv) + ((u.ind.juvenile + u.ind.adult) * p.ind.K_X_juv)
    u.ind.f_X = (u.glb.X / p.glb.V_patch) / ((u.glb.X / p.glb.V_patch) + K_X)

    dI_emb = u.ind.embryo * (Complex(u.ind.S + u.ind.E_mt)^(2/3)).re * p.ind.dI_max_emb * u.ind.y_T
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
function maintenance_M2!(du, u, p, t)::Nothing 

    k_M = (u.embryo + u.larva + u.metamorph) * p.k_M_emb + (u.juvenile + u.adult) * p.k_M_juv
    k_J = (u.embryo + u.larva + u.metamorph) * p.k_J_emb + (u.juvenile + u.adult) * p.k_J_juv

    du.M = (u.S + u.E_mt) * k_M * u.y_j[2] * u.y_jP[2] * u.y_T # somatic maintenance flux
    du.J = u.H * k_J * u.y_j[2] * u.y_jP[2] * u.y_T # maturity maintenance flux

    return nothing
end


"""
    growth_M2!(du, u, p, t)::Tuple{Real,Real}

Calculation of growth fluxes for amphibians. <br>
Returns life stage-specific values for `eta_AS` and `kappa`.
"""
function growth_M2!(du, u, p, t, eta_AS, kappa)::Nothing
    
    #### somatic growth for embryos, juveniles and adults ####
    # apply shrinking equation if maintenance costs are not covered 
    # we use a sigmoid switch to avoid yet another if/else statement

    dS_emb_juv_ad = sig( 
        kappa * du.A, 
        du.M, 
        -(du.M / p.eta_SA - kappa * du.A), # shrinking equation applies if kappa*dA < dM
        u.y_j[1] * eta_AS * (kappa * du.A - du.M) # otherwise, "normal" growth (but what is normal, really?)
    )    
    
    ##### somatic growth for larvae ####
    # growth for larvae follows the same rule as for juveniles and adults, except that there is an additional split in the resource allocation
    # as before, shrinking equation 
    
    dS_lrv = u.larva * sig( 
        kappa * du.A, 
        du.M, 
        -(du.M / p.eta_SA - (1 - p.gamma) * kappa * du.A), 
        eta_AS * u.y_j[1] * ((1 - p.gamma) * kappa * du.A - du.M)
    )

    #### somatic growth for metamorphs ####
    # for metamorph, structural growth is assumed to be driven by the residual assimilation flux

    dS_mt = u.metamorph * eta_AS * u.y_j[1] * u.y_jP[1] * du.A  
    du.S = (u.embryo + u.juvenile + u.adult) * dS_emb_juv_ad + u.larva * dS_lrv + u.metamorph * dS_mt

    return nothing
end 


"""
    metamorphic_reserve!(du, u, p, t, kappa)::Nothing

Calculation of metamorphic reserve dynamics for amphibians. <br>
Metamorphic reserve is accumulated during larval development and depleted during metamorphic climax. 
"""
function metamorphic_reserve_M2!(du, u, p, t, eta_AS, kappa)::Nothing
    
    # the metamorphic reserve is fueled by the gamma*kappa-fraction of the assimilation flux
    dE_mt_lrv = eta_AS * u.y_j[1] * (p.gamma * kappa * du.A - du.M) # before metamorphosis, reserve is built up
    dE_mt_mt = -(du.H + du.J + du.M) # during metamorphosis, while there is still E_mt left, it will be used to fuel maintenance and maturation
    du.E_mt = u.larva * dE_mt_lrv + u.metamorph * dE_mt_mt
    du.E_mt_max = u.larva * du.E_mt

    return nothing
end

""" 
    Amphibian_DEB_M2!(du, u, p, t)::Nothing

Complete ODE system for alternative formulation of the AmphiDEB model (see `Amphibian_DEB_M2!`)

"""
function AmphiDEB_ODE_M2!(du, u, p, t)::Nothing

    EcotoxSystems.DEBODE_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    AmphiDEB_individual_M2!(du, u, p, t)

    return nothing
end


function AmphiDEB_individual_M2!(du, u, p, t)::Nothing

    TKTD_mix_IA!(du, u, p, t) # TKTD following default model
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    Amphibian_DEB_M2!(du, u, p, t) # Amphibian DEB model

    return nothing
end

""" 
    Amphibian_DEB_M2!(du, u, p, t)::Nothing

Alternative formulation of the AmphiDEB model, where `E_mt` is viewed as a sub-compartment of structure.
This leads to the following main differences:

- `E_mt` is subject to somatic maintenance 
- `E_mt` contributes to the surface area scaling factor
- `E_mt` is affected by growth efficiency, as well as effects on growth efficiency

"""
function Amphibian_DEB_M2!(du, u, p, t)::Nothing

    determine_life_stage!(du.ind, u.ind, p.ind, t)
    y_T!(du, u, p, t)
    eta_AS, kappa = life_stage_effects(du, u, p, t)

    ingestion_M2!(du, u, p, t)
    maintenance_M2!(du.ind, u.ind, p.ind, t)
    growth_M2!(du.ind, u.ind, p.ind, t, eta_AS, kappa)
    maturation!(du.ind, u.ind, p.ind, t, kappa)
    metamorphic_reserve_M2!(du.ind, u.ind, p.ind, t, eta_AS, kappa)

    reproduction!(du.ind, u.ind, p.ind, t, kappa)

    return nothing
end
