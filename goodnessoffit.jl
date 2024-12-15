using StatsBase, Random
   
function r_mudi(K::Int, N, mu, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    mc_max_guess = Int(round(N*mu*fit_m*eff*10))
    p_draws = p_mudi(mc_max_guess, N, mu, fit_m, eff)
    k = 1
    cumulative_p = p_draws[k]
    uni_draws = sort(rand(K))
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p
            k += 1
            if k > mc_max_guess
                mc_max_guess *= 10
                p_draws = p_mudi(mc_max_guess, N, mu, fit_m, eff)
            end
            cumulative_p += p_draws[k]
        end
        random_draws[j] = k-1
    end
    return random_draws, p_draws[1:k]
end

function r_mudi(K::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    mc_max_guess = Int(round(N*mu_off*((1-f_on)*fit_m + S*f_on)*eff*10))
    p_draws = p_mudi(mc_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    k = 1
    cumulative_p = p_draws[k]
    uni_draws = sort(rand(K))
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p
            k += 1
            if k > mc_max_guess
                mc_max_guess *= 10
                p_draws = p_mudi(mc_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
            end
            cumulative_p += p_draws[k]
        end
        random_draws[j] = k-1
    end
    return random_draws, p_draws[1:k]
end

function Hellinger_dist(R::Int, num_c::Int, mc_draws, p_draws)
    H = Vector{Float64}(undef, R)
    s = randperm(R * num_c) 
    c_indices = reshape(s, num_c, R)
    for i in eachindex(H)
        mc_sample = mc_draws[c_indices[:, i]]
        p_sample = counts(mc_sample, 0:maximum(mc_sample)) ./ num_c
        p = p_draws[1:maximum(mc_sample)+1]
        norm_factor = sqrt(sum(p_sample) * sum(p))
        H[i] = sqrt(1 - sum(sqrt.(p_sample .* p)) / norm_factor)
    end
    return H
end

function gof(mc_counts, p_counts, R::Int, num_c, mc_draws, p_draws)
    H = Hellinger_dist(R, num_c, mc_draws, p_draws)
    p_sample = mc_counts ./ num_c
    norm_factor = sqrt(sum(p_sample) * sum(p_counts))
    H_sample = sqrt(1 - sum(sqrt.(p_sample .* p_counts)) / norm_factor)
    return H_sample, ecdf(H)(H_sample)
end

    