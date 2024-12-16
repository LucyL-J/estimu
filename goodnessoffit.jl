using StatsBase

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

    