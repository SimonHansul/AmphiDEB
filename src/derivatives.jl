# derivatives.jl
# model functions which are used by all models

## alternative model configurations


"""
    AmphiDEB_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

ODE system with log-logistic toxicodynamics. 
"""
function AmphiDEB_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    AmphiDEB_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    AmphiDEB_individual_ODE_with_loglogistic_TD!(du, u, p, t)

    return nothing
end

# linear toxicodynamics (classic "DEBtox" stress function)

"""
    AmphiDEB_ODE_with_linear_TD!(du, u, p, t)::Nothing

ODE system with liner toxicodynamics. 
"""
function AmphiDEB_ODE_with_linear_TD!(du, u, p, t)::Nothing

    AmphiDEB_global!(du, u, p, t)
    Pathogen_growth!(du, u, p, t)
    AmphiDEB_individual_ODE_with_linear_TD!(du, u, p, t)

    return nothing
end

### individual-level ODE configurations

function AmphiDEB_individual_ODE_with_loglogistic_TD!(du, u, p, t)::Nothing

    TKTD_LL2!(du, u, p, t) # TKTD with mixtures assuming IA
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    Amphibian_DEB!(du, u, p, t) # Amphibian DEB model

    return nothing
end

function AmphiDEB_individual_ODE_with_linear_TD!(du, u, p, t)::Nothing

    TKTD_linear!(du, u, p, t) # TKTD with mixtures assuming IA
    Pathogen_Infection!(du, u, p, t) # infection, release of zoospores and relative response to sporangia density
    Amphibian_DEB!(du, u, p, t) # Amphibian DEB mode

    return nothing

end

""" 
    Amphibian_DEB!(du, u, p, t)::Nothing

The default amphibian DEB model assumes a metamorphic reserve compartment,
which is accumulated during larval development, and depleted during metamorphosis. 
During metamorphosis, ingestion rate decreases gradually and reaches 0 at the end of metamorphosis.
"""
function Amphibian_DEB!(du, u, p, t)::Nothing

    determine_life_stage!(du, u, p, t)
    Arrhenius!(du, u, p, t)
    eta_AS, kappa = life_stage_and_plasticity_effects(du, u, p, t)

    ingestion!(du, u, p, t)
    maintenance!(du, u, p, t)
    growth!(du, u, p, t, eta_AS, kappa)
    maturation!(du, u, p, t, kappa)
    metamorphic_reserve!(du, u, p, t, eta_AS, kappa)

    reproduction!(du, u, p, t, kappa)

    return nothing
end

## callbacks (events)

# adding pathogen spores at a fixed simulation time
condition_inoculation(u, t, integrator) = integrator.p.glb.pathogen_inoculation_time - t
function effect_inoculation!(integrator) 
    integrator.u.glb.P_Z = integrator.p.glb.pathogen_inoculation_dose
end

# medium renewal, result in removal of spores
condition_renewal(u, t, integrator) = prod(integrator.p.glb.medium_renewals  .- t)  
function effect_renewal!(integrator)
    integrator.u.glb.P_Z = 0.
end

# function to generate define all callbacks in a single CallbackSet
function AmphODE_callbacks()

    cb_inoculation = ContinuousCallback(condition_inoculation, effect_inoculation!)
    cb_renewal = ContinuousCallback(condition_renewal, effect_renewal!)

    return CallbackSet(
        cb_inoculation, 
        cb_renewal
        )

end

## derivative functions

### global derivatives

function AmphiDEB_global!(du, u, p, t)::Nothing

    for i in eachindex(u.glb[:X])
        du.glb.X[i] = dX(p.glb[:dX_in][i], p.glb[:k_V][i], u.glb[:X][i])
    end

    #u.ind.T = p.glb.tempfun(t, p.glb.temp...) # TODO: temporally variable temperature could be added here, e.g. applying temperature_sinusoidal

    return nothing
end


"""
    temperature_sinusoidal(t::Float64, T_max::Float64, T_min::Float64, t_peak ::Float64)::Float64

Calculate seasonal fluctuations in temperature from a sinusoidal function. 
"""
@inline function temperature_sinusoidal(t::Float64, T_mean::Float64, T_amp::Float64, T_phi::Float64)::Float64

    return T_mean + T_amp * sin*(2π /365 * t + T_phi)

end

@inline function dX(
    dX_in::Float64, 
    k_V::Float64, 
    X::Float64
    )::Float64

    return dX_in - k_V * X

end

"""
    function dP_Z(
        mu::Float64,
        P_Z::Float64
        )::Float64

Derivative of the zoospore abundance `P_Z`. 

This function does not take changes in zoopsore abundance due to infection into account. 
This is handled in `Pathogen_Infection!`.

args: 

- `mu`: Zoospore background mortality rate
- `P`_Z`: Current zoospore abundance in the environment
"""
@inline function dP_Z(
    mu::Float64,
    P_Z::Float64
    )::Float64

    return -mu * P_Z

end

function Pathogen_growth!(du, u, p, t)::Nothing
    
    du.glb.P_Z = dP_Z(p.pth[:mu], u.glb[:P_Z]) # change in zoospores in environment
    
    return nothing
end

## individual-level derivatives

"""
    function dP_S(
        v0::Float64,
        gamma::Float64,
        P_Z::Float64,
        eta::Float64,
        f::Float64,
        P_S::Float64,
        sigma0::Float64,
        sigma1::Float64,
        Chi::Float64
        )::Float64

Derivative of the individual-specific sporangia abundace `P_S`.
"""
@inline function dP_S(
    v0::Float64,
    gamma::Float64,
    P_Z::Float64,
    eta::Float64,
    f::Float64,
    P_S::Float64,
    sigma0::Float64,
    sigma1::Float64,
    Chi::Float64
    )::Float64

    return v0 * gamma * P_Z + v0 * eta * f * P_S - (sigma0 + sigma1 * Chi * P_S) * P_S

end

function Pathogen_Infection!(du, u, p, t)::Nothing

    # life-stage specificity of the zoospore-host encounter rate 
    # - currently turned off to simplify things 
    # this needs to be re-evaluated. tadpoles do get infected but are unlikely to die 
    # maybe only make the effect parameters life stage-specific?
    gamma = p.pth.gamma 

    # growth and killing of sporangia 
    du.ind.P_S = dP_S(p.pth[:v0], p.pth[:gamma], u.glb[:P_Z], p.pth[:eta], p.pth[:f], u.ind[:P_S], p.pth[:sigma0], p.pth[:sigma1], p.ind[:Chi])

    # feedback with zoospore population 
    du.glb.P_Z += p.pth.eta * (1-p.pth.f) * u.ind.P_S # release of spores
    du.glb.P_Z -= gamma * u.glb.P_Z # encystment

    # relative response to pathogen 
   
    @. u.ind.y_jP = LL2(u.ind.P_S, p.ind.E_P, p.ind.B_P) 
    u.ind.y_jP[2] /= u.ind.y_jP[2]^2 # converting a monotonically decreasing to increasing response

    return nothing
end

#### TKTD functions

@inline function LL2(x::Float64, p::NTuple{2,Float64})::Float64
    return (1 / (1 + Complex(x / p[1]) ^ p[2])).re
end

@inline function LL2pos(x::Float64, p::NTuple{2,Float64})::Float64
    return 1 - log((1 / (1 + Complex(x / p[1]) ^ p[2])).re)
end

@inline function NEC2pos(x::Float64, p::NTuple{2,Float64})::Float64
    return 1 + (p[2] * max(0, x - p[1]))
end

@inline function NEC2neg(x::Float64, p::NTuple{2,Float64})::Float64
    return min(1, 1 - (p[2] * max(0, x - p[1])))
end

@inline LL2(x::Float64, p1::Float64, p2::Float64)::Float64 = LL2(x, (p1, p2))
@inline LL2pos(x::Float64, p1::Float64, p2::Float64)::Float64 = LL2pos(x, (p1, p2))
@inline LL2GUTS(x::Float64, p1::Float64, p2::Float64)::Float64 = -log(LL2(x, (p1, p2)))


@inline NEC2neg(x::Float64, p1::Float64, p2::Float64)::Float64 = NEC2neg(x, (p1, p2))
@inline NEC2pos(x::Float64, p1::Float64, p2::Float64)::Float64 = NEC2pos(x, (p1, p2))


@inline function minimal_TK(
    embryo::Float64,
    k_D::Float64, 
    C_W::Float64,
    D::Float64
    )::Float64

    return (1-embryo) * k_D * (C_W - D)

end

"""
    TKTD_LL2!(du, u, p, t)::Nothing

TKTD model with following configuration: 

- Mixture toxicity based on independent action (IA)
- Log-logistic relationship between damage and metabolic processes
"""
@inline function TKTD_LL2!(du, u, p, t)::Nothing
    
    @unpack glb, ind = u

    ind.y_j .= 1.0 # reset relative responses 
    ind.h_z = p.ind[:h_b] # reset GUTS-SD hazard rate to background mortality

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(ind.y_j) # for every PMoA
            # calculate change in damage
            du.ind.D_j[z,j] = minimal_TK(ind.embryo, p.ind.KD[z,j], glb.C_W[z], ind.D_j[z,j]) 
            # update relative response with respect to PMoA j
            # PMoAs with decreasing response
            if !(y in [2,6]) 
                ind.y_j[j] *= LL2(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            # PMoAs with increasing response
            else
                ind.y_j[j] *= LL2pos(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.ind.D_h[z] = minimal_TK(ind[:embryo], p.ind[:KD_h][z], glb[:C_W][z], ind[:D_h][z]) 
        # update hazard rate
        ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

"""
    TKTD_linear!(du, u, p, t)::Nothing

TKTD model with following configuration: 

- Mixture toxicity based on independent action (IA)
- Linear-above-threshold relationship between damage and metabolic processes (the default DEBtox stress function)
"""
@inline function TKTD_linear!(du, u, p, t)::Nothing

    @unpack glb, ind = u

    ind.y_j .= 1.0 # reset relative responses 
    ind.h_z = p.ind[:h_b] # reset GUTS-SD hazard rate to background mortality

    for z in eachindex(glb.C_W) # for every chemical
        for j in eachindex(ind.y_j) # for every PMoA
            # calculate change in damage
            du.ind.D_j[z,j] = minimal_TK(ind.embryo, p.ind.KD[z,j], glb.C_W[z], ind.D_j[z,j]) 
            # update relative response with respect to PMoA j
            # PMoAs with decreasing response
            if j != 2 
                ind.y_j[j] *= NEC2neg(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            # PMoAs with increasing response
            else
                ind.y_j[j] *= NEC2pos(ind.D_j[z,j], p.ind.E[z,j], p.ind.B[z,j])
            end
        end
        # calculate change in damage for lethal effects
        du.ind.D_h[z] = minimal_TK(ind[:embryo], p.ind[:KD_h][z], glb[:C_W][z], ind[:D_h][z]) 
        # update hazard rate
        ind.h_z += LL2GUTS(ind.D_h[z], p.ind.E_h[z], p.ind.B_h[z])
    end

    du.ind.S_z = -ind.h_z * ind.S_z # survival probability according to GUTS-RED-SD
    
    return nothing
end

@inline function is_embryo(
    X_emb::Float64;
    beta::Float64 = 1e3
    )::Float64

    return sig(X_emb, 0., 0., 1., beta = beta)

end

@inline function is_larva(
    X_emb::Float64,
    H::Float64,
    H_j1::Float64,
    y_H_neg::Float64, 
    y_H_pos::Float64;
    beta1::Float64 = 1e3,
    beta2::Float64 = 1e7
    )::Float64

    return sig(X_emb, 0., 1., 0., beta = beta1) * sig(H, H_j1 * y_H_neg * y_H_pos, 1., 0., beta = beta2)

end

@inline function is_metamorph(
    H::Float64, 
    H_j1::Float64, 
    y_H_neg::Float64,
    y_H_pos::Float64,
    E_mt::Float64;
    beta1::Float64 = 1e7, 
    beta2::Float64 = 1e3, 
    )::Float64

    return sig(H, H_j1 * y_H_neg * y_H_pos, 0., 1., beta = beta1) * sig(E_mt, 0., 0., 1., beta = beta2) 

end

@inline function is_juvenile(
    H::Float64,
    H_j1::Float64,
    y_H_neg::Float64,
    y_H_pos::Float64,
    E_mt::Float64,
    H_p::Float64;
    beta1 = 1e3,
    beta2 = 1e3,
    beta3 = 1e3
    )::Float64

    return sig(H, H_j1 * y_H_neg * y_H_pos, 0., 1., beta = beta1) * sig(E_mt, 0., 1., 0., beta = beta2) * sig(H, H_p, 1., 0., beta = beta3)

end

@inline function is_adult(
    H::Float64,
    H_p::Float64;
    beta = 1e3
    )::Float64

    return sig(H, H_p, 0., 1., beta = beta)

end


function determine_life_stage!(du, u, p, t)::Nothing

    u.ind.embryo = is_embryo(u.ind[:X_emb])
    u.ind.larva = is_larva(u.ind[:X_emb], u.ind[:H], p.ind[:H_j1], u.ind[:y_j][5], u.ind[:y_j][6])
    u.ind.metamorph = is_metamorph(u.ind[:H], p.ind[:H_j1], u.ind[:y_j][5], u.ind[:y_j][6], u.ind[:E_mt]) 
    u.ind.juvenile = is_juvenile(u.ind[:H], p.ind[:H_j1], u.ind[:y_j][5], u.ind[:y_j][6], u.ind[:E_mt], p.ind[:H_p]) 
    u.ind.adult = is_adult(u.ind[:H], p.ind[:H_p]) 

    return nothing
end

@inline function calc_eta_AS(
    embryo::Float64,
    larva::Float64,
    metamorph::Float64,
    eta_AS_emb::Float64,
    juvenile::Float64,
    adult::Float64,
    eta_AS_juv::Float64
    )::Float64

    return (embryo + larva + metamorph) * eta_AS_emb + (juvenile + adult) * eta_AS_juv

end


"""
    calc_kappa(
        embryo::Float64,
        larva::Float64,
        metamorph::Float64,
        kappa_emb::Float64,
        juvenile::Float64,
        adult::Float64,
        kappa_juv::Float64,
        b_T::Float64, 
        T_ref::Float64,
        y_K::Float64
        )::Float64

Calculate `kappa` accounting for life-stage, temperature effects and chemical effects. 

Temperature effects are implemented according to Romoli et al. (2024).
"""
@inline function calc_kappa(
    embryo::Float64,
    larva::Float64,
    metamorph::Float64,
    kappa_emb::Float64,
    juvenile::Float64,
    adult::Float64,
    kappa_juv::Float64,
    b_T::Float64, 
    T_ref::Float64,
    T::Float64,
    y_K::Float64
    )::Float64

    kappa = (embryo + larva + metamorph) * kappa_emb + (juvenile + adult) * kappa_juv

    return 1/(1 + (((1-kappa)/kappa) * exp(-b_T * ((T_ref - T)/T_ref))))*y_K
end


"""
    life_stage_and_plasticity_effects(du, u, p, t)::Tuple{Float64,Float64}

Handles life-stage specificity of parameters. 
The parameter notation foresees that the superscript indicates the first life stage for which a value is valid,
e.g. if we have eta_AS_emb and eta_AS_juv, then eta_AS_emb is applied for eymbros, larvae and metamorphs, and eta_AS_juv is applied for juveniles and adults
"""
function life_stage_and_plasticity_effects(du, u, p, t)::Tuple{Float64,Float64}

    eta_AS = calc_eta_AS(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:eta_AS_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:eta_AS_juv])
    kappa = calc_kappa(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:kappa_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:kappa_juv], p.ind[:b_T], p.ind[:T_ref], p.glb[:T], u.ind.y_j[7])

    return eta_AS, kappa
end

"""
    y_T(
        T_A::Float64,
        T_ref::Float64,
        T::Float64
        )::Float64

Calculates temperature correction coefficient `y_T` according to Arrhenius equation.

- `T_A` : Arrhenius temperature (K)
- `T_ref` : Reference temperature (K)
- `T` : Current ambient temperature
"""
@inline function y_T(
    T_A::Float64,
    T_ref::Float64,
    T::Float64
    )::Float64

    return exp((T_A / T_ref) - (T_A / T))

end

function Arrhenius!(du, u, p, t)::Nothing
    
    u.ind.y_T = y_T(p.ind[:T_A], p.ind[:T_ref], p.glb[:T])

    return nothing

end


@inline function f_X(
    X::Float64,
    V_patch::Float64,
    K_X::Float64
    )::Float64

    return (X / V_patch) / ((X / V_patch) + K_X)

end


@inline function f_X(
    larva::Float64,
    metamorph::Float64, 
    juvenile::Float64, 
    adult::Float64,
    X::Vector{Float64},
    V_patch::Union{Vector{Real},Vector{Float64}}, 
    K_X_lrv::Float64,
    K_X_juv::Float64
    )::Float64


    return ((larva + metamorph) * f_X(X[1], V_patch[1], K_X_lrv)) + ((juvenile + adult) * f_X(X[2], V_patch[2], K_X_juv))

end

@inline function calc_dI_emb(
    embryo::Float64,
    S::Float64,
    dI_max_emb::Float64,
    y_T::Float64
    )::Float64

    return embryo * (Complex(S)^(2/3)).re * dI_max_emb * y_T
    
end

@inline function calc_dI_mt(
    metamorph::Float64,
    f_X::Float64,
    dI_max_lrv::Float64,
    E_mt::Float64,
    E_mt_max::Float64,
    S::Float64,
    y_T::Float64
    )::Float64

    return metamorph * f_X * dI_max_lrv * (E_mt / E_mt_max) * (Complex(S)^(2/3)).re * y_T

end

@inline function calc_dI_lrv(
    larva::Float64,
    f_X::Float64,
    dI_max_lrv::Float64,
    S::Float64,
    y_T::Float64,
    )::Float64

    return larva * f_X * dI_max_lrv * (Complex(S)^(2/3)).re * y_T

end

@inline function calc_dI_juv(
    juvenile::Float64,
    adult::Float64,
    f_X::Float64,
    dI_max_juv::Float64,
    S::Float64,
    y_T::Float64
    )::Float64

    return (juvenile + adult) * f_X * dI_max_juv * (Complex(S)^(2/3)).re * y_T

end

@inline function dI(
    dI_emb::Float64,
    dI_mt::Float64, 
    dI_lrv::Float64, 
    dI_juv_ad::Float64
    )::Float64

    return dI_emb + dI_mt + dI_lrv + dI_juv_ad

end

@inline function dA(
    dI::Float64, 
    eta_IA::Float64, 
    y_A::Float64, 
    y_AP::Float64
    )::Float64

    return dI * eta_IA * y_A * y_AP

end

"""
    ingestion!(du, u, p, t)::Nothing

Life stage-specific calculation of ingestion rate for amphibians. <br>
Note that this function requires all components to be included in the arguments `du`, `u` and `p`, so that we can access the external food concentration. <br>
Hence, we have `u.ind`, `u.glb`, `p.ind`, etc., instead of simply `u` and `p`.
"""
function ingestion!(du, u, p, t)::Nothing

    u.ind.f_X = f_X(u.ind[:larva], u.ind[:metamorph], u.ind[:juvenile], u.ind[:adult], u.glb[:X], p.glb[:V_patch], p.ind[:K_X_lrv], p.ind[:K_X_juv])

    dI_emb = calc_dI_emb(u.ind[:embryo], u.ind[:S], p.ind[:dI_max_emb], u.ind[:y_T]) # ingestion by embryos
    dI_mt = calc_dI_mt(u.ind[:metamorph], u.ind[:f_X], p.ind[:dI_max_lrv], u.ind[:E_mt], u.ind[:E_mt_max], u.ind[:S], u.ind[:y_T]) # residual ingestion flux by metamorphs (which may include tadpoles close to climax)
    dI_lrv = calc_dI_lrv(u.ind[:larva], u.ind[:f_X], p.ind[:dI_max_lrv], u.ind[:S], u.ind[:y_T]) # ingestion by larvae
    dI_juv_ad = calc_dI_juv(u.ind[:juvenile], u.ind[:adult], u.ind[:f_X], p.ind[:dI_max_juv], u.ind[:S], u.ind[:y_T])
    
    du.ind.I = dI(dI_emb, dI_mt, dI_lrv, dI_juv_ad)
    du.ind.X_emb = -dI_emb
    du.glb.X[1] -= (u.ind[:larva] + u.ind[:metamorph]) * du.ind.I
    du.glb.X[2] -= (u.ind[:juvenile] + u.ind[:adult]) * du.ind.I

    # assimilation flux
    du.ind.A = dA(du.ind[:I], p.ind[:eta_IA], u.ind[:y_j][3], u.ind[:y_jP][3])

    return nothing 
end

@inline function determine_k_M(
    embryo::Float64,
    larva::Float64,
    metamorph::Float64,
    k_M_emb::Float64,
    juvenile::Float64,
    adult::Float64,
    k_M_juv::Float64
    )::Float64

    return (embryo + larva + metamorph) * k_M_emb + (juvenile + adult) * k_M_juv

end

@inline function determine_k_J(
    embryo::Float64,
    larva::Float64,
    metamorph::Float64,
    k_J_emb::Float64,
    juvenile::Float64,
    adult::Float64,
    k_J_juv::Float64
    )::Float64

    return (embryo + larva + metamorph) * k_J_emb + (juvenile + adult) * k_J_juv

end

@inline function dM(
    S::Float64,
    k_M::Float64,
    y_M::Float64,
    y_MP::Float64,
    y_T::Float64
    )::Float64

    return S * k_M * y_M * y_MP * y_T

end

@inline function dJ(
    H::Float64,
    k_J::Float64,
    y_M::Float64,
    y_MP::Float64,
    y_T::Float64
    )::Float64

    return H * k_J * y_M * y_MP * y_T

end

"""
    maintenance!(du, u, p, t)::Nothing 

Calculation of maintenance fluxes for amphibians.
"""
function maintenance!(du, u, p, t)::Nothing 

    k_M = determine_k_M(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:k_M_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:k_M_juv])   
    k_J = determine_k_J(u.ind[:embryo], u.ind[:larva], u.ind[:metamorph], p.ind[:k_J_emb], u.ind[:juvenile], u.ind[:adult], p.ind[:k_J_juv])   

    du.ind.M = dM(u.ind[:S], k_M, u.ind[:y_j][2], u.ind[:y_jP][2], u.ind[:y_T])
    du.ind.J = dJ(u.ind[:H], k_J, u.ind[:y_j][2], u.ind[:y_jP][2], u.ind[:y_T])

    return nothing
end

@inline function dR(
    adult::Float64,
    eta_AR::Float64,
    y_R::Float64,
    y_RP::Float64,
    kappa::Float64,
    dA::Float64,
    dJ::Float64
    )::Float64

    return adult * clipneg(eta_AR * y_R * y_RP * ((1 - kappa) * dA - dJ))

end

function reproduction!(du, u, p, t, kappa::Float64)::Nothing
    
    du.ind.R = dR(u.ind[:adult], p.ind[:eta_AR] , u.ind[:y_j][4], u.ind[:y_jP][4], kappa, du.ind[:A], du.ind[:J])

    return nothing
end

@inline function calc_dS_emb_juv_ad(
    kappa::Float64,
    dA::Float64,
    dM::Float64,
    eta_SA::Float64,
    y_G::Float64,
    y_GP::Float64,
    eta_AS::Float64
    )::Float64

    return sig( 
        kappa * dA, 
        dM, 
        -(dM / eta_SA - kappa * dA), 
        y_G * y_GP * eta_AS * (kappa * dA - dM) 
    )    

end

@inline function calc_dS_lrv(
    kappa::Float64,
    dA::Float64,
    dM::Float64,
    eta_SA::Float64,
    gamma::Float64,
    eta_AS::Float64,
    y_G::Float64,
    y_GP::Float64
    )::Float64

    if (kappa * dA) > dM
        return eta_AS * y_G * y_GP * (1 - gamma) * (kappa * dA - dM)
    else
        return -(dM / eta_SA - (1 - gamma) * kappa * dA)
    end
end

@inline function calc_dS_mt(
    metamorph::Float64,
    eta_AS::Float64,
    y_G::Float64,
    y_GP::Float64,
    dA::Float64
    )::Float64

    return metamorph * eta_AS * y_G * y_GP * dA

end

@inline function dS(
    embryo::Float64,
    juvenile::Float64,
    adult::Float64,
    dS_emb_juv_ad::Float64,
    larva::Float64,
    dS_lrv::Float64,
    metamorph::Float64,
    dS_mt::Float64
    )::Float64

    return (embryo + juvenile + adult) * dS_emb_juv_ad + larva * dS_lrv + metamorph * dS_mt

end

"""
    growth!(du, u, p, t)

Calculation of life stage-specific growth fluxes for amphibians. <br>
"""
function growth!(du, u, p, t, eta_AS::Float64, kappa::Float64)::Nothing
    
    dS_emb_juv_ad = calc_dS_emb_juv_ad(kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], u.ind[:y_j][1], u.ind[:y_jP][1], eta_AS)
    dS_lrv = calc_dS_lrv(kappa, du.ind[:A], du.ind[:M], p.ind[:eta_SA], p.ind[:gamma], eta_AS, u.ind[:y_j][1], u.ind[:y_jP][1])
    dS_mt = calc_dS_mt(u.ind[:metamorph], eta_AS, u.ind[:y_j][1], u.ind[:y_jP][1], du.ind[:A])
    
    du.ind.S = dS(u.ind[:embryo], u.ind[:juvenile], u.ind[:adult], dS_emb_juv_ad, u.ind[:larva], dS_lrv, u.ind[:metamorph], dS_mt)

    return nothing
end 

# metamorphic reserve dynamics for larvae
@inline function calc_dE_mt_lrv(
    eta_AS::Float64, 
    y_G::Float64,
    y_G_P::Float64,
    gamma::Float64,
    kappa::Float64,
    dA::Float64,
    dM::Float64
    )::Float64

    return eta_AS * y_G * y_G_P * gamma * (kappa * dA - dM)

end

# metamorphic reserve dynamics for metamorphs
@inline function calc_dE_mt_mt(
    dH::Float64,
    dJ::Float64,
    dM::Float64
    )::Float64

    return -(dH + dJ + dM)

end

# metamorphic reserve dynamics for any life stage
@inline function dE_mt(
    larva::Float64,
    dE_mt_lrv::Float64,
    metamorph::Float64,
    dE_mt_mt::Float64
    )::Float64

    return larva * dE_mt_lrv + metamorph * dE_mt_mt

end

# tracking maximum reserve level
@inline function dE_mt_max(
    larva::Float64,
    dE_mt::Float64
    )::Float64

    return larva * dE_mt

end

"""
    metamorphic_reserve!(du, u, p, t, kappa)::Nothing

Calculation of metamorphic reserve dynamics for amphibians. <br>
Metamorphic reserve is accumulated during larval development and depleted during metamorphic climax. 
"""
function metamorphic_reserve!(du, u, p, t, eta_AS::Float64, kappa::Float64)::Nothing
    
    dE_mt_lrv = calc_dE_mt_lrv(p.ind[:eta_AS_emb], u.ind[:y_j][1], u.ind[:y_jP][1], p.ind[:gamma], kappa, du.ind[:A], du.ind[:M]) 
    dE_mt_mt = calc_dE_mt_mt(du.ind[:H], du.ind[:J], du.ind[:M])
    du.ind.E_mt = dE_mt(u.ind[:larva], dE_mt_lrv, u.ind[:metamorph], dE_mt_mt)
    du.ind.E_mt_max = dE_mt_max(u.ind[:larva], du.ind[:E_mt])

    return nothing
end

@inline function calc_dH(
    kappa::Float64,
    dA::Float64,
    dJ::Float64
    )::Float64

    return clipneg((1 - kappa) * dA - dJ)

end

@inline function calc_dH_mt(
    kappa::Float64,
    dM::Float64,
    dJ::Float64
    )::Float64

    return clipneg((1 - kappa) * dM) / kappa - dJ

end

@inline function dH(
    adult::Float64,
    metamorph::Float64,
    dH::Float64,
    dH_mt::Float64
    )::Float64

    return (1 - adult) * (1 - metamorph) * dH + metamorph * dH_mt
end

"""
    maturation!(du, u, p, t, kappa)::Nothing

Calculation of maturation fluxes for amphibians.
"""
function maturation!(du, u, p, t, kappa::Float64)::Nothing

    # maturation follows kappa-rule for all but metamorphs
    dH_all = calc_dH(kappa, du.ind[:A], du.ind[:J])

    ## for metamorhps, we assume that all but somatic growth is fueled by E_mt, 
    # allowing us to apply the kappa rule to calculate dH from dM and dJ
    dH_mt = calc_dH_mt(kappa, du.ind[:M], du.ind[:J])

    # for adults, there is no maturation => we only need to differentiate between non-adults and metamorph
    du.ind.H = dH(u.ind[:adult], u.ind[:metamorph], dH_all, dH_mt)

    return nothing
end
