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
    u_1 = find_zero(LL_ratio, (m, 10*mc_max))
    return [l_1 u_1]
end
function CI_m_fitm(mc_counts, mc_max, m, inv_fit_m, eff, MLL)
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
    u_1 = find_zero(LL_ratio_1, (m, 10*mc_max))
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
    u = find_zero(LL_ratio, (m, 10*mc_max))
    return [l u]
end
# Homogeneous response with jointly inferred mutant fitness
function CI_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff, MLL)
    M = Vector{Float64}(undef, 6)
    function LL_ratio_1(para)
        if para == m_UT
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para, P[1], P[2], eff)
            res = Optim.optimize(log_likelihood_1, [m_S, inv_fit_m])
            P_res = Optim.minimizer(res)
            return -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para, P_res[1], P_res[2], eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m_UT))
    l_1_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, l_1, P[1], P[2], eff)
    res = Optim.optimize(l_1_M, [m_S, inv_fit_m])
    M[1] = Optim.minimizer(res)[1]/l_1
    u_1 = find_zero(LL_ratio_1, (m_UT, 10*mc_max_UT))
    u_1_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, u_1, P[1], P[2], eff)
    res = Optim.optimize(u_1_M, [m_S, inv_fit_m])
    M[2] = Optim.minimizer(res)[1]/u_1
    function LL_ratio_2(para)
        if para == m_S
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], para, P[2], eff)
            res = Optim.optimize(log_likelihood_2, [m_UT, inv_fit_m])
            P_res = Optim.minimizer(res)
            return -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P_res[1], para, P_res[2], eff) - MLL - chisq_1_95/2
        end
    end
    l_2 = find_zero(LL_ratio_2, (0., m_S))
    l_2_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], l_2, P[2], eff)
    res = Optim.optimize(l_2_M, [m_UT, inv_fit_m])
    M[3] = l_2/Optim.minimizer(res)[1]
    u_2 = find_zero(LL_ratio_2, (m_S, 10*mc_max_S))
    u_2_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], u_2, P[2], eff)
    res = Optim.optimize(u_2_M, [m_UT, inv_fit_m])
    M[4] = u_2/Optim.minimizer(res)[1]
    function LL_ratio_3(para)
        if para == inv_fit_m
            return -chisq_1_95/2
        else
            log_likelihood_3(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], P[2], para, eff)
            res = Optim.optimize(log_likelihood_3, [m_UT, m_S])
            P_res = Optim.minimizer(res)
            return -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P_res[1], P_res[2], para, eff) - MLL - chisq_1_95/2
        end
    end
    l_3 = find_zero(LL_ratio_3, (0., inv_fit_m))
    l_3_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], P[2], l_3, eff)
    res = Optim.optimize(l_3_M, [m_UT, m_S])
    M[5] = Optim.minimizer(res)[2]/Optim.minimizer(res)[1]
    u_3 = Inf
    try
        u_3 = find_zero(LL_ratio_3, (inv_fit_m, Inf))
    catch err
    end
    if u_3 == Inf
        eff = [e[1] for e in eff]
        u_3_M_0(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, P[1], P[2])
        res = Optim.optimize(u_3_M_0, [m_UT, m_S].*eff)
        Optim.minimizer(res) .*= eff
    else
        u_3_M(P) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], P[2], u_3, eff)
        res = Optim.optimize(u_3_M, [m_UT, m_S])
    end
    M[6] = Optim.minimizer(res)[2]/Optim.minimizer(res)[1]
    return [l_1 u_1; l_2 u_2; l_3 u_3; minimum(M) maximum(M)]
end