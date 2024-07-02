function log_likelihood_m(mc_counts, mc_max, m, q0, q)
    if m <= 0.
        return -Inf
    else
        p = mudi(mc_max, m, q0, q)
        @views ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end