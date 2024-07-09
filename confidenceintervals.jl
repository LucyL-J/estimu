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
function CI_m_fitm(mc_counts, mc_max, m, inv_fit_m, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_m_fitm(mc_counts, mc_max, para, P[1])
            res = Optim.optimize(log_likelihood_1, [inv_fit_m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, para, P1) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(LL_ratio_1, (m, mc_max))
    catch err
        u_1 = find_zero(LL_ratio_1, (m, 10*mc_max))
    end
    function LL_ratio_2(para) 
        if para == inv_fit_m
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_m_fitm(mc_counts, mc_max, P[1], para)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, P1, para) - MLL - chisq_1_95/2
        end
    end
    l_2 = find_zero(LL_ratio_2, (0., inv_fit_m))
    u_2 = Inf
    try
        u_2 = find_zero(LL_ratio_2, (inv_fit_m, Inf))
    catch err
    end
    return [l_1 u_1; l_2 u_2]
end
function CI_m_fitm(mc_counts, mc_max, m, inv_fit_m, MLL, eff)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_m_fitm(mc_counts, mc_max, para, P[1], eff)
            res = Optim.optimize(log_likelihood_1, [inv_fit_m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, para, P1, eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(LL_ratio_1, (m, mc_max))
    catch err
        u_1 = find_zero(LL_ratio_1, (m, 10*mc_max))
    end
    function LL_ratio_2(para) 
        if para == inv_fit_m
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_m_fitm(mc_counts, mc_max, P[1], para, eff)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, P1, para, eff) - MLL - chisq_1_95/2
        end
    end
    l_2 = find_zero(LL_ratio_2, (0., inv_fit_m))
    u_2 = Inf
    try
        u_2 = find_zero(LL_ratio_2, (inv_fit_m, Inf))
    catch err
    end
    return [l_1 u_1; l_2 u_2]
end
function CI_m_fitm(mc_counts, mc_max, m, inv_fit_m, MLL, eff, small_eff::Bool)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_m_fitm(mc_counts, mc_max, para, P[1], eff, true)
            res = Optim.optimize(log_likelihood_1, [inv_fit_m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, para, P1, eff, true) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = m
    try
        u_1 = find_zero(LL_ratio_1, (m, mc_max))
    catch err
        u_1 = find_zero(LL_ratio_1, (m, 10*mc_max))
    end
    function LL_ratio_2(para) 
        if para == inv_fit_m
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_m_fitm(mc_counts, mc_max, P[1], para, eff, true)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_m_fitm(mc_counts, mc_max, P1, para, eff, true) - MLL - chisq_1_95/2
        end
    end
    l_2 = find_zero(LL_ratio_2, (0., inv_fit_m))
    u_2 = Inf
    try
        u_2 = find_zero(LL_ratio_2, (inv_fit_m, Inf))
    catch err
    end
    return [l_1 u_1; l_2 u_2]
end
# Without change in mutation rate (fixed mutant fitness, single or multiple stressful cond.)
function CI_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, q0_UT, q_UT, q0_S, q_S, MLL)
    function LL_ratio(para)
        if para == m
            return -chisq_1_95/2
        else
            return -log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, q0_UT, q_UT, q0_S, q_S) - MLL - chisq_1_95/2
        end
    end
    l = find_zero(LL_ratio, (0., m))
    u = m
    try
        u = find_zero(LL_ratio, (m, mc_max))
    catch err
        u = find_zero(LL_ratio, (m, 10*mc_max))
    end
    return [l u]
end