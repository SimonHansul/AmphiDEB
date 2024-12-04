function Amphibian_DEB_structural_mobilization!(du, u, p, t)::Nothing
    """ 
    The structural mobilization model assumes that metamorphs use structure to fuel 
    maturation and maintenance. Ingestion stops immediately as individuals
    """

    @unpack glb,ind,pth = u 

    # temperature correction
    ind.y_T = exp((p.ind.T_A / p.ind.T_ref) - (p.ind.T_A / p.glb.T)) 

    # ingestion rates and feedback with resource pools
    K_X = (u.ind.larva * u.ind.metamorph * p.ind.K_X_lrv) + (u.ind.juvenile * u.ind.adult * p.ind.K_X_juv)
    ind.f_X = (glb.X / p.glb.V_patch) / ((glb.X / p.glb.V_patch) + K_X)
    
    #### life stage-specific resource uptake ####
    
    dI_emb = ind.embryo * (Complex(ind.S)^(2/3)).re * p.ind.dI_max_emb * ind.y_T
    dI_mt = 0 #ind.metamorph * ind.f_X * p.ind.dI_max_lrv * (u.ind.E_mt / u.ind.E_mt_max) * (Complex(ind.S)^(2/3)).re * ind.y_T
    dI_lrv = ind.larva * ind.f_X * p.ind.dI_max_lrv * (Complex(ind.S)^(2/3)).re * ind.y_T
    dI_juv_ad = (ind.juvenile + ind.adult) * ind.f_X * p.ind.dI_max_juv * (Complex(ind.S)^(2/3)).re * ind.y_T
    
    du.ind.I = dI_emb + dI_mt + dI_lrv + dI_juv_ad
    du.ind.X_emb = -dI_emb

    # assimilation flux
    du.ind.A = du.ind.I * p.ind.eta_IA * u.ind.y_j[3] * u.ind.y_jP[3]

    # maintenance fluxes
    k_M = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.k_M_emb + (u.ind.juvenile + u.ind.adult) * p.ind.k_M_juv
    k_J = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.k_J_emb + (u.ind.juvenile + u.ind.adult) * p.ind.k_J_juv

    du.ind.M = ind.S * k_M * ind.y_j[2] * ind.y_jP[2] * ind.y_T # Somatic maintenance flux
    du.ind.J = ind.H * k_J * ind.y_j[2] * ind.y_jP[2] * ind.y_T # Maturity maintenance flux

    #### life-stage specific growth flux ###

    # Somatic growth for embryos, juveniles and adults follows the default model
    eta_AS = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.eta_AS_emb + (u.ind.juvenile + u.ind.adult) * p.ind.eta_AS_juv
    kappa = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.kappa_emb + (u.ind.juvenile + u.ind.adult) * p.ind.kappa_juv

    # maturation
    dH = clipneg(((1 - kappa) * du.ind.A) - du.ind.J) 
    dH_mt = ((1 - kappa) * du.ind.M) / kappa - du.ind.J
    du.ind.H = (1 - u.ind.adult) * dH + u.ind.metamorph * dH_mt

    # somatic growth
    dS = sig( 
        kappa * du.ind.A, 
        du.ind.M, -(du.ind.M / p.ind.eta_SA - kappa * du.ind.A), 
        ind.y_j[1] * eta_AS * (kappa * du.ind.A - du.ind.M)
    )   

    dS_mt = -(du.M + du.H + du.J)/p.ind.eta_SA
    du.ind.S = (1 - ind.metamorph) * dS + ind.metamorph * dS_mt


   
    # reproduction flux
    du.ind.R = ind.adult * clipneg(p.ind.eta_AR * ind.y_j[4] * ind.y_jP[4] * ((1 - kappa) * du.ind.A - du.ind.J))  # reproduction flux

    return nothing
end

function Amphibian_DEB_soft_structural_mobilization!(du, u, p, t)::Nothing
    """ 
    In the soft structural mobilization model, 
    structure is used during metamorphosis, 
    but ingestion decreases gradually (depending on maturity). 
    The residual ingestion flux buffers the decrease in structure.
    """

    @unpack glb,ind,pth = u 

    # temperature correction
    ind.y_T = exp((p.ind.T_A / p.ind.T_ref) - (p.ind.T_A / p.glb.T)) 

    # ingestion rates and feedback with resource pools
    K_X = (u.ind.larva * u.ind.metamorph * p.ind.K_X_lrv) + (u.ind.juvenile * u.ind.adult * p.ind.K_X_juv)
    ind.f_X = (glb.X / p.glb.V_patch) / ((glb.X / p.glb.V_patch) + K_X)
    
    #### life stage-specific resource uptake ####
    
    dI_emb = ind.embryo * (Complex(ind.S)^(2/3)).re * p.ind.dI_max_emb * ind.y_T
    dI_mt = ind.metamorph * ind.f_X * p.ind.dI_max_lrv * ((p.ind.H_46 * u.ind.y_j[5] * u.ind.y_jP[5])/ind.H) * (Complex(ind.S)^(2/3)).re * ind.y_T
    dI_lrv = ind.larva * ind.f_X * p.ind.dI_max_lrv * (Complex(ind.S)^(2/3)).re * ind.y_T
    dI_juv_ad = (ind.juvenile + ind.adult) * ind.f_X * p.ind.dI_max_juv * (Complex(ind.S)^(2/3)).re * ind.y_T
    
    du.ind.I = dI_emb + dI_mt + dI_lrv + dI_juv_ad
    du.ind.X_emb = -dI_emb

    # assimilation flux
    du.ind.A = du.ind.I * p.ind.eta_IA * u.ind.y_j[3] * u.ind.y_jP[3]

    # maintenance fluxes
    k_M = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.k_M_emb + (u.ind.juvenile + u.ind.adult) * p.ind.k_M_juv
    k_J = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.k_J_emb + (u.ind.juvenile + u.ind.adult) * p.ind.k_J_juv

    du.ind.M = ind.S * k_M * ind.y_j[2] * ind.y_jP[2] * ind.y_T # Somatic maintenance flux
    du.ind.J = ind.H * k_J * ind.y_j[2] * ind.y_jP[2] * ind.y_T # Maturity maintenance flux

    #### life-stage specific growth flux ###

    # Somatic growth for embryos, juveniles and adults follows the default model
    eta_AS = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.eta_AS_emb + (u.ind.juvenile + u.ind.adult) * p.ind.eta_AS_juv
    kappa = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.kappa_emb + (u.ind.juvenile + u.ind.adult) * p.ind.kappa_juv

    # maturation
    dH = clipneg(((1 - kappa) * du.ind.A) - du.ind.J) 
    dH_mt = ((1 - kappa) * du.ind.M) / kappa - du.ind.J
    du.ind.H = (1 - u.ind.adult) * dH + u.ind.metamorph * dH_mt

    # somatic growth
    dS = sig( 
        kappa * du.ind.A, 
        du.ind.M, -(du.ind.M / p.ind.eta_SA - kappa * du.ind.A), 
        ind.y_j[1] * eta_AS * (kappa * du.ind.A - du.ind.M)
    )   

    dS_mt = -((du.M + du.H + du.J)/p.ind.eta_SA - du.ind.A)
    du.ind.S = (1 - ind.metamorph) * dS + ind.metamorph * dS_mt

    # reproduction flux
    du.ind.R = ind.adult * clipneg(p.ind.eta_AR * ind.y_j[4] * ind.y_jP[4] * ((1 - kappa) * du.ind.A - du.ind.J))  # reproduction flux

    return nothing
end

# a maturity threshold is used to determine the end of metamorphosis
condition_juvenile_structural_mobilization(u, t, integrator) = u.ind.H - integrator.p.ind.H_46 * u.ind.y_j[5] * u.ind.y_jP[5]

function structural_mobilization_callbacks()

    cb_larva = ContinuousCallback(condition_larva, effect_larva!)
    cb_metamorph = ContinuousCallback(condition_metamorph, effect_metamorph!)
    cb_juvenile = ContinuousCallback(condition_juvenile_structural_mobilization, effect_juvenile!)
    cb_adult = ContinuousCallback(condition_adult, effect_adult!)

    cb_inoculation = ContinuousCallback(condition_inoculation, effect_inoculation!)
    cb_renewal = ContinuousCallback(condition_renewal, effect_renewal!)

    return CallbackSet(
        cb_larva, 
        cb_metamorph, 
        cb_juvenile, 
        cb_adult, 
        cb_inoculation, 
        cb_renewal
        )
end

