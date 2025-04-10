function robustsum(v::AbstractVector)::Float64

    if length(v)>0
        return sum(v)
    else
        return 0
    end

end


"""
    AmphiDEB_global_rules!(m)

Global rule-based portion of the AmphiDEB IBM. 
"""
function AmphiDEB_global_rules!(m)
    m.u.glb.N = length(m.individuals) # tracking population size
    m.u.glb.X = max.(0, m.u.glb.X) # HOTFIX : negative resource abundances can cause chaos
    
    # record counts per life stage
    m.u.glb.N_emb = robustsum(map(x -> isapprox(1, x.u.ind.embryo, atol = 0.01), m.individuals))
    m.u.glb.N_lrv = robustsum(map(x -> isapprox(1, x.u.ind.larva, atol = 0.01), m.individuals))
    m.u.glb.N_mt = robustsum(map(x -> isapprox(1, x.u.ind.metamorph, atol = 0.01), m.individuals))
    m.u.glb.N_juv = robustsum(map(x -> isapprox(1, x.u.ind.juvenile, atol = 0.01), m.individuals))
    m.u.glb.N_ad = robustsum(map(x -> isapprox(1, x.u.ind.adult, atol = 0.01), m.individuals))
    
end
