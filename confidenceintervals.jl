using Roots

# Profile log-likelihood to calculate confidence intervals
# The 95% quantile of the Chi-squared distribution used to calculate CIs
# Depend on the observed mutant counts, on which parameters are estimated, the ML estimates and the maximum likelihood itself 
chisq_1_95 = 3.84145882069412447634704221854917705059051513671875
# Stadard model (with optional diff. mutant fitness)
function CI_m(mc_counts, mc_max, m, q0, q, MLL)
    function LL_ratio(para) 
        if para == m
            return -chisq_1_95/2
        else
            return -log_likelihood_m(mc_counts, mc_max, para, q0, q) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio, (0., m))
    u_1 = m
    try
        u_1 = find_zero(LL_ratio, (m, mc_max))
    catch err
        u_1 = find_zero(LL_ratio, (m, 10*mc_max))
    end
    return [l_1 u_1]
end