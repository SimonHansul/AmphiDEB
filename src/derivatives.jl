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
   
    @. u.ind.y_jP = 1. # EcotoxSystems.LL2(u.ind.P_S/(Complex(u.ind.S^(2/3)).re), p.ind.e_P, p.ind.b_P) 
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

# mixture TKTD based on independent action model
@inline function TKTD_mix_IA!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    # scaled damage dynamics based on the minimal model

    #@. du.ind.D_z = (1 - ind.embryo) * p.ind.k_D_z * (glb.C_W - ind.D_z)
    #@. du.ind.D_h = (1 - ind.embryo) * p.ind.k_D_h * (glb.C_W - ind.D_h)

    for z in eachindex(glb.C_W)
        # for sublethal effects, we broadcost over all PMoAs
        @. du.ind.D_z[z,:] = (1 - ind.embryo) * @view(p.ind.k_D_z[z,:]) * (glb.C_W[z] - @view(ind.D_z[z,:]))
        # for lethal effects, we have only one value per stressor
        du.ind.D_h[z] = (1 - ind.embryo) * p.ind.k_D_h[z] * (glb.C_W[z] - ind.D_h[z])
    end

    @. ind.y_z = EcotoxSystems.softNEC2neg(ind.D_z, p.ind.e_z, p.ind.b_z) # relative responses per stressor and PMoA
    
    ind.y_j .= reduce(*, ind.y_z; dims=1) # relative responses per PMoA are obtained as the product over all chemical stressors
    ind.y_j[2] /= ind.y_j[2]^2 # for pmoas with increasing responses (M), the relative response has to be inverted  (x/x^2 == 1/x) 

    #ind.h_z = sum(@. softNEC2GUTS(ind.D_h, p.ind.e_h, p.ind.b_h)) # hazard rate according to GUTS-RED-SD
    ind.h_z = 0 
    @inbounds for z in eachindex(ind.D_h)
        ind.h_z += EcotoxSystems.softNEC2GUTS(ind.D_h[z], p.ind.e_h[z], p.ind.b_h[z])
    end
    
    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

# definition of life stages with sigmoid functions
function determine_life_stage!(du, u, p, t)::Nothing

    u.embryo = sig(u.X_emb, 0, 0, 1, beta = 1e3) # if there is still vitellus left, we are in embryo stage
    u.larva = sig(u.X_emb, 0, 1, 0, beta = 1e3) * sig(u.H, p.H_j1 * u.y_j[5], 1, 0, beta = 1e7) # if the embryo is used up but the next maturity threshold is not reached, we are in larval stage
    u.metamorph =  sig(u.H, p.H_j1 * u.y_j[5], 0, 1, beta = 1e7) * sig(u.E_mt, 0, 0, 1, beta = 1e3) # above the maturity threshold for Gosner stage 42, while there is still metamorphic reserve left, we are in metamorph stage
    u.juvenile = sig(u.H, p.H_j1 * u.y_j[5], 0, 1, beta = 1e3) * sig(u.E_mt, 0, 1, 0, beta = 1e3) * sig(u.H, p.H_p, 1, 0, beta = 1e3) # after metamorphosis but below the threshold for puberty, we are in juvenile stage
    u.adult = sig(u.H, p.H_p, 0, 1, beta = 1e3) # adult stage is reached beyond maturity threshold H_p

    return nothing
end


"""
    life_stage_effects(du, u, p, t)::Tuple{Float64,Float64}

Handles life-stage specificity of parameters. 
The parameter notation foresees that the superscript indicates the first life stage for which a value is valid,
e.g. if we have eta_AS_emb and eta_AS_juv, then eta_AS_emb is applied for eymbros, larvae and metamorphs, and eta_AS_juv is applied for juveniles and adults
"""
function life_stage_effects(du, u, p, t)::Tuple{Float64,Float64}

    eta_AS = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.eta_AS_emb + (u.ind.juvenile + u.ind.adult) * p.ind.eta_AS_juv
    kappa = (u.ind.embryo + u.ind.larva + u.ind.metamorph) * p.ind.kappa_emb + (u.ind.juvenile + u.ind.adult) * p.ind.kappa_juv

    # adjust kappa for temperature
    kappa = 1/(1 + (((1-kappa)/kappa) * exp(-p.ind.b_T * ((p.ind.T_ref - p.glb.T)/p.ind.T_ref))))

    return eta_AS,kappa
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

"""
    growth!(du, u, p, t)::Tuple{Real,Real}

Calculation of growth fluxes for amphibians. <br>
Returns life stage-specific values for `eta_AS` and `kappa`.
"""
function growth!(du, u, p, t, eta_AS, kappa)::Nothing
    
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
        eta_AS * u.y_j[1] * (1 - p.gamma) * ( kappa * du.A - du.M)
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
function metamorphic_reserve!(du, u, p, t, eta_AS, kappa)::Nothing
    
    # the metamorphic reserve is fueled by the gamma*kappa-fraction of the assimilation flux
    dE_mt_lrv = p.gamma * (kappa * du.A - du.M) # before metamorphosis, reserve is built up
    dE_mt_mt = -(du.H + du.J + du.M) # during metamorphosis, while there is still E_mt left, it will be used to fuel maintenance and maturation
    du.E_mt = u.larva * dE_mt_lrv + u.metamorph * dE_mt_mt
    du.E_mt_max = u.larva * du.E_mt

    return nothing
end

"""
    maturation!(du, u, p, t, kappa)::Nothing

Calculation of maturation fluxes for amphibians.
"""
function maturation!(du, u, p, t, kappa)::Nothing

    # maturation follows kappa-rule for all but metamorphs
    dH = clipneg((1 - kappa) * du.A - du.J)

    ## for metamorhps, we assume that all but somatic growth is fueled by E_mt, 
    # allowing us to apply the kappa rule to calculate dH from dM and dJ
    dH_mt = clipneg((1 - kappa) * du.M) / kappa - du.J

    # for adults, there is no maturation => we only need to differentiate between non-adults and metamorph
    du.H = (1 - u.adult) * (1 - u.metamorph) * dH + u.metamorph * dH_mt

    return nothing
end



function AmphiDEB_ODE!(du, u, p, t)::Nothing

    EcotoxSystems.DEBODE_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    AmphiDEB_individual!(du, u, p, t)

    return nothing
end

@inline function temperature_sinusoidal(t::Float64, T_max::Float64, T_min::Float64, t_peak ::Float64)::Float64

    amplitude = (T_max - T_min) / 2
    offset = (T_max + T_min) / 2
    omega = 2π / 365
    phase_shift = (π / 2) - omega * t_peak
    return amplitude * sin(omega * t + phase_shift) + offset

end

function AmphiDEB_global!(du, u, p, t)::Nothing

    EcotoxSystems.DEBODE_global!(du, u, p, t)
    #u.glb.T = p.glb.tempfun(t, p.glb.temp...) 

    return nothing
end


function AmphiDEB_individual!(du, u, p, t)::Nothing

    TKTD_mix_IA!(du, u, p, t) # TKTD following default model
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    Amphibian_DEB!(du, u, p, t) # Amphibian DEB model

    return nothing
end

""" 
    Amphibian_DEB!(du, u, p, t)::Nothing

The default amphibian DEB model assumes a metamorphic reserve compartment,
which is accumulated during larval development, and depleted during metamorphosis. 
During metamorphosis, ingestion rate decreases gradually and reaches 0 at the end of metamorphosis.
"""
function Amphibian_DEB!(du, u, p, t)::Nothing

    determine_life_stage!(du.ind, u.ind, p.ind, t)
    y_T!(du, u, p, t)
    eta_AS, kappa = life_stage_effects(du, u, p, t)

    ingestion!(du, u, p, t)
    maintenance!(du.ind, u.ind, p.ind, t)
    growth!(du.ind, u.ind, p.ind, t, eta_AS, kappa)
    maturation!(du.ind, u.ind, p.ind, t, kappa)
    metamorphic_reserve!(du.ind, u.ind, p.ind, t, eta_AS, kappa)

    reproduction!(du.ind, u.ind, p.ind, t, kappa)

    return nothing
end
