using StatsBase

function LL_dist(R::Int, num_c::Int, N, mu, fit_m, eff)
    mc_rdraws, p_rdraws = r_mudi(R*num_c, N, mu, fit_m, eff)
    LL = Vector{Float64}(undef, R)
    s = randperm(R * num_c) 
    c_indices = reshape(s, num_c, R)
    st = false
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        if any(x -> x == 0, p)
            println("Warning: Zero probability in p_rdraws at index $i")
            println("p: $p")
        end
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
        if isnan(LL[i])
            LL[i] = Inf  # Assign a large value to avoid NaN
            st = true
        end
    end
    if any(st)
        println("Warning: NaN values found in LL; mu = $mu, fit_m = $fit_m")
    end
    return LL
end

function LL_dist(R::Int, num_c::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    mc_rdraws, p_rdraws = r_mudi(R*num_c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    LL = Vector{Float64}(undef, R)
    s = randperm(R * num_c) 
    c_indices = reshape(s, num_c, R)
    st = false
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        if any(x -> x == 0, p)
            println("Warning: Zero probability in p_rdraws at index $i")
            println("p: $p")
        end
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
        if isnan(LL[i])
            LL[i] = Inf  # Assign a large value to avoid NaN
            st = true
        end
    end
    if any(st)
        println("Warning: NaN values found in LL; mu_off = $mu_off, S = $S")
    end
    return LL
end