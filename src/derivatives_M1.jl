# derivatives_M1.jl 
# the M1 model assumes that the metamorphic reserve fuels obligatory processes during metamorphic climax (maturation + maintenance fluxes), 
# the residual ingestion flux fuels facultative processes (somatic growth)
# model equations for model version AmphiDEB M1

function AmphiDEB_ODE_M1!(du, u, p, t)::Nothing

    EcotoxSystems.DEBODE_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    AmphiDEB_individual_M1!(du, u, p, t)

    return nothing
end

function AmphiDEB_individual_M1!(du, u, p, t)::Nothing

    EcotoxSystems.TKTD_mix_IA!(du, u, p, t) # TKTD following default model
    Pathogen_Infection!(du, u, p, t) 
    Amphibian_DEB_M1!(du, u, p, t)

    return nothing
end

"""
    growth!(du, u, p, t)::Tuple{Real,Real}

Calculation of growth fluxes for amphibians. <br>
Returns life stage-specific values for `eta_AS` and `kappa`.
"""
function growth_M1!(du, u, p, t)::Tuple{Real,Real}
    
    #### handle life-stage specificity of parameters ####
    # the parameter notation foresees that the superscript indicates the first life stage for which a value is valid
    # e.g. if we have eta_AS_emb and eta_AS_juv, then eta_AS_emb is applied for eymbros, larvae and metamorphs, and eta_AS_juv is applied for juveniles and adults

    eta_AS = (u.embryo + u.larva + u.metamorph) * p.eta_AS_emb + (u.juvenile + u.adult) * p.eta_AS_juv
    kappa = (u.embryo + u.larva + u.metamorph) * p.kappa_emb + (u.juvenile + u.adult) * p.kappa_juv

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

    return eta_AS, kappa
end 

"""
    maturation!(du, u, p, t, kappa)::Nothing

Calculation of maturation fluxes for amphibians.
"""
function maturation_M1!(du, u, p, t, kappa)::Nothing

    # maturation follows kappa-rule for all but metamorphs
    dH = clipneg((1 - kappa) * du.A - du.J)

    ## for metamorhps, we assume that all but somatic growth is fueled by E_mt, 
    # allowing us to apply the kappa rule to calculate dH from dM and dJ
    dH_mt = clipneg((1 - kappa) * du.M) / kappa - du.J

    # for adults, there is no maturation => we only need to differentiate between non-adults and metamorph
    du.H = (1 - u.adult) * (1 - u.metamorph) * dH + u.metamorph * dH_mt

    return nothing
end

"""
    metamorphic_reserve!(du, u, p, t, kappa)::Nothing

Calculation of metamorphic reserve dynamics for amphibians. <br>
Note that the metamorphic reserve does not work the same as reserve in the standard DEB model. <br>
Metamorphic reserve is accumulated during larval development and depleted during metamorphic climax. 
"""
function metamorphic_reserve_M1!(du, u, p, t, eta_AS, kappa)::Nothing
    
    # the metamorphic reserve is fueled by the gamma*kappa-fraction of the assimilation flux
    dE_mt_lrv = (1 - p.gamma) * kappa * du.A - du.M # before metamorphosis, reserve is built up
    dE_mt_mt = -(du.H + du.J + du.M) # during metamorphosis, while there is still E_mt left, it will be used to fuel maintenance and maturation
    du.E_mt = u.larva * dE_mt_lrv + u.metamorph * dE_mt_mt
    du.E_mt_max = u.larva * du.E_mt

    return nothing
end


""" 
    Amphibian_DEB!(du, u, p, t)::Nothing

The default amphibian DEB model assumes a metamorphic reserve compartment,
which is accumulated during larval development, and depleted during metamorphosis. 
During metamorphosis, ingestion rate decreases gradually and reaches 0 at the end of metamorphosis.
"""
function Amphibian_DEB_M1!(du, u, p, t)::Nothing

    determine_life_stage!(du.ind, u.ind, p.ind, t)
    y_T!(du, u, p, t)

    ingestion!(du, u, p, t)
    maintenance!(du.ind, u.ind, p.ind, t)
    eta_AS, kappa = growth_M1!(du.ind, u.ind, p.ind, t)
    maturation_M1!(du.ind, u.ind, p.ind, t, kappa)
    metamorphic_reserve_M1!(du.ind, u.ind, p.ind, t, eta_AS, kappa)

    reproduction!(du.ind, u.ind, p.ind, t, kappa)

    return nothing
end
