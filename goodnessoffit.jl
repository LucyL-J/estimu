using StatsBase, Random
   
function r_mudi(K::Int, N, mu, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    mc_max_guess = Int(round(N*mu*fit_m*eff*10))
    probs = p_mudi(mc_max_guess, N, mu, fit_m, eff)
    k = 1
    cumulative_p = probs[k]
    uni_draws = sort(rand(K))
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p
            k += 1
            if k > mc_max_guess
                mc_max_guess *= 10
                probs = p_mudi(mc_max_guess, N, mu, fit_m, eff)
            end
            cumulative_p += probs[k]
        end
        random_draws[j] = k-1
    end
    return random_draws
end

function r_mudi(K::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    mc_max_guess = Int(round(N*mu_off*((1-f_on)*fit_m + S*f_on)*eff*10))
    probs = p_mudi(mc_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    k = 1
    cumulative_p = probs[k]
    uni_draws = sort(rand(K))
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p
            k += 1
            if k > mc_max_guess
                mc_max_guess *= 10
                probs = p_mudi(mc_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
            end
            cumulative_p += probs[k]
        end
        random_draws[j] = k-1
    end
    return random_draws
end

#mc_total = r_mudi(R*c, N, mu, fit_m, eff)
#probs = p_mudi(maximum(mc_total), N, mu, fit_m, eff)
#mc_total = r_mudi(R*c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
#probs = p_mudi(maximum(mc_total), N, mu_off, S, f_on, rel_div_on, fit_m, eff)
function Hellinger_dist(R, c, mc_total, probs)
    H = Vector{Float64}(undef, R)
    s = randperm(R * c) 
    c_indices = reshape(s, c, R)
    for i in eachindex(H)
        mc_sample = mc_total[c_indices[:, i]]
        p_sample = counts(mc_sample, 0:maximum(mc_sample)) ./ c
        p = probs[collect(1:maximum(mc_sample)+1)]
        norm_factor = sqrt(sum(p_sample) * sum(p))
        H[i] = sqrt(1 - sum(sqrt.(p_sample .* p)) / norm_factor)
    end
    return H
end

function gof_measure(mc_sample, p, H)
    p_sample = counts(mc_sample, 0:maximum(mc_sample)) ./ length(mc_sample)
    norm_factor = sqrt(sum(p_sample) * sum(p))
    H_sample = sqrt(1 - sum(sqrt.(p_sample .* p)) / norm_factor)
    return H_sample, ecdf(H)(H_sample)
end

    