# check final structural mass
function default_individual_rules!(
    a::EcotoxSystems.AbstractDEBIndividual, 
    m::EcotoxSystems.AbstractDEBIBM
    )::Nothing

    @unpack glb,ind = a.u
    du = a.du
    p = a.p

    ind.age += m.dt

    # death due to aging
    if ind.age > p.ind.a_max
        ind.cause_of_death = 1.
    end

    # life-stage transitions are part of ODE in amphibian model and omitted here
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    ind.S_max_hist = max(ind.S, ind.S_max_hist)

    if ((ind.S/ind.S_max_hist) < p.ind.S_rel_crit) && (rand() > exp(-p.ind.h_S * m.dt))
        ind.cause_of_death = 2.
    end

    # reproduction, assuming a constant reproduction period

    # reproduction only occurs if the reproduction period has been exceeded
    if ind.time_since_last_repro >= p.ind.tau_R
        # if that is the case, calculate the number of offspring,
        # based on the reproduction buffer and the dry mass of an egg
        let num_offspring = trunc(ind.R / p.ind.X_emb_int)
            for _ in 1:num_offspring
                m.idcount += 1 # increment individual counter
                push!(m.individuals, EcotoxSystems.DEBIndividual( # create new individual and push to individuals vector
                    m.p,
                    m.u.glb;
                    id = m.idcount,
                    cohort = Int(ind.cohort + 1),
                    individual_ode! = a.individual_ode!,
                    individual_rules! = a.individual_rules!,
                    init_individual_statevars = a.init_individual_statevars,
                    gen_ind_params = a.generate_individual_params
                    )
                )
                ind.R -= p.ind.X_emb_int # decrease reproduction buffer
                ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
            end
            ind.time_since_last_repro = 0. # reset reproduction period
        end
    # if reproduction period has not been exceeded,
    else
        ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end