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
function log_likelihood_m_fitm(mc_counts, mc_max, m, inv_fit_m)
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = -1
        q = q_coeffs(mc_max, inv_fit_m)
        p = mudi(mc_max, m, q0, q)
        @views ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_fitm(mc_counts, mc_max, m, inv_fit_m, eff)
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = q0_coeff(inv_fit_m, eff)
        q = q_coeffs(mc_max, inv_fit_m, eff)
        p = mudi(mc_max, m, q0, q)
        @views ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_fitm(mc_counts, mc_max, m, inv_fit_m, eff, small_eff::Bool)
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = q0_coeff(inv_fit_m, eff)
        q = q_coeffs(mc_max, inv_fit_m, eff, true)
        p = mudi(mc_max, m, q0, q)
        @views ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end