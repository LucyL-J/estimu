using StatsBase

function LL_dist(num_c::Int, N, mu, fit_m, eff)
    mc_rdraws, p_rdraws = r_mudi(R_gof*num_c, N, mu, fit_m, eff)
    LL = Vector{Float64}(undef, R_gof)
    s = randperm(R_gof * num_c) 
    c_indices = reshape(s, num_c, R_gof)
    st = false
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
    end
    return LL, length(p_rdraws), p_rdraws[end]
end

function LL_dist(num_c::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    mc_rdraws, p_rdraws = r_mudi(R_gof*num_c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    LL = Vector{Float64}(undef, R_gof)
    s = randperm(R_gof * num_c) 
    c_indices = reshape(s, num_c, R_gof)
    st = false
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
    end
    return LL, length(p_rdraws), p_rdraws[end]
end