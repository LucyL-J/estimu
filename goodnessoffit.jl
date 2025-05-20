using StatsBase

function Hellinger_dist(R::Int, num_c::Int, mc_rdraws, p_rdraws)
    H = Vector{Float64}(undef, R)
    s = randperm(R * num_c) 
    c_indices = reshape(s, num_c, R)
    for i in eachindex(H)
        mc_observed = mc_rdraws[c_indices[:, i]]
        p_observed = counts(mc_observed, 0:maximum(mc_observed)) ./ num_c
        p = p_rdraws[1:maximum(mc_observed)+1]
        norm_factor = sqrt(sum(p_observed) * sum(p))
        H[i] = sqrt(1 - sum(sqrt.(p_observed .* p)) / norm_factor)
    end
    return H
end

function gof(mc_counts, p_counts, R::Int, num_c, mc_rdraws, p_rdraws)
    H = Hellinger_dist(R, num_c, mc_rdraws, p_rdraws)
    p_observed = mc_counts ./ num_c
    norm_factor = sqrt(sum(p_observed) * sum(p_counts))
    H_observed = sqrt(1 - sum(sqrt.(p_observed .* p_counts)) / norm_factor)
    return H_observed, ecdf(H)(H_observed)
end

    