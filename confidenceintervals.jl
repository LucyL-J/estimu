using Roots

# Profile log-likelihood to calculate confidence intervals
# The 95% quantile of the Chi-squared distribution used to calculate CIs
# Depend on the observed mutant counts, on which parameters are estimated, the ML estimates and the maximum likelihood itself 
chisq_1_95 = 3.84145882069412447634704221854917705059051513671875
# Standard model (fixed mutant fitness)
function CI_m(mc_counts, mc_max, m, q0, q, eff, MLL)
    function LL_ratio(para) 
        if para == m
            return -chisq_1_95/2
        else
            return -log_likelihood_m(mc_counts, mc_max, para, q0, q) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio, (0., m))
    u_1 = find_zero(LL_ratio, (m, m_max(mc_max,eff)))
    return [l_1 u_1]
end
# Standard model (mutant fitness inferred)
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
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max,eff)))
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
# Without change in mutation rate (fixed mutant fitness)
function CI_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, q0_UT, q_UT, q0_S, q_S, eff, MLL)
    function LL_ratio(para)
        if para == m
            return -chisq_1_95/2
        else
            return -log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, q0_UT, q_UT, q0_S, q_S) - MLL - chisq_1_95/2
        end
    end
    l = find_zero(LL_ratio, (0., m))
    u = find_zero(LL_ratio, (m, m_max(mc_max,eff)))
    return [l u]
end
# Without change in mutation rate (jointly inferred mutant fitness)
function CI_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, inv_fit_m, eff, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para, P[1], eff)
            res = Optim.optimize(log_likelihood_1, [inv_fit_m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para, P1, eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max,eff)))
    function LL_ratio_2(para)
        if para == inv_fit_m
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P[1], para, eff)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P1, para, eff) - MLL - chisq_1_95/2
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
# Homogeneous response (fixed mutant fitness)
function CI_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, q0_UT, q_UT, q0_S, q_S, MLL)
    function min_M(P)
        if -log_likelihood_m(mc_counts_UT, mc_max_UT, P[1], q0_UT, q_UT) - log_likelihood_m(mc_counts_S, mc_max_S, P[2], q0_S, q_S) - MLL > chisq_1_95/2
            return Inf
        else
            return P[2]/P[1]
        end
    end
    res = Optim.optimize(min_M, [m_UT, m_S])
    l = Optim.minimum(res)
    function max_M(P)
        if -log_likelihood_m(mc_counts_UT, mc_max_UT, P[1], q0_UT, q_UT) - log_likelihood_m(mc_counts_S, mc_max_S, P[2], q0_S, q_S) - MLL > chisq_1_95/2
            return Inf
        else
            return -P[2]/P[1]
        end
    end
    res = Optim.optimize(max_M, [m_UT, m_S])
    u = -Optim.minimum(res)
    return [l u]
end
# Homogeneous response (mutant fitness inferred separately)
function CI_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, inv_fit_m_UT, inv_fit_m_S, eff, MLL)
    eff_UT, eff_S = eff
    if eff[1] == 1.
        eff_UT = false
    elseif eff[1] < 0.5
        eff_UT = (eff[1], true)
    end
    if eff[2] == 1.
        eff_S = false
    elseif eff[2] < 0.5
        eff_S = (eff[2], true)
    end
    eff = (eff_UT, eff_S)
    function min_M(P)
        if -log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, P[1], P[3], eff[1]) - log_likelihood_m_fitm(mc_counts_S, mc_max_S, P[2], P[4], eff[2]) - MLL > chisq_1_95/2
            return Inf
        else
            return P[2]/P[1]
        end
    end
    res = Optim.optimize(min_M, [m_UT, m_S, inv_fit_m_UT, inv_fit_m_S])
    l_1 = Optim.minimum(res)
    function max_M(P)
        if -log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, P[1], P[3], eff[1]) - log_likelihood_m_fitm(mc_counts_S, mc_max_S, P[2], P[4], eff[2]) - MLL > chisq_1_95/2
            return Inf
        else
            return -P[2]/P[1]
        end
    end
    res = Optim.optimize(max_M, [m_UT, m_S, inv_fit_m_UT, inv_fit_m_S])
    u_1 = -Optim.minimum(res)
    function min_ifit_ratio(P)
        if -log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, P[1], P[3], eff[1]) - log_likelihood_m_fitm(mc_counts_S, mc_max_S, P[2], P[4], eff[2]) - MLL > chisq_1_95/2
            return Inf
        else
            return P[4]/P[3]
        end
    end
    res = Optim.optimize(min_ifit_ratio, [m_UT, m_S, inv_fit_m_UT, inv_fit_m_S])
    l_2 = Optim.minimum(res)
    function max_ifit_ratio(P)
        if -log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, P[1], P[3], eff[1]) - log_likelihood_m_fitm(mc_counts_S, mc_max_S, P[2], P[4], eff[2]) - MLL > chisq_1_95/2
            return Inf
        else
            return -P[4]/P[3]
        end
    end
    res = Optim.optimize(max_ifit_ratio, [m_UT, m_S, inv_fit_m_UT, inv_fit_m_S])
    u_2 = -Optim.minimum(res)
    return [l_1 u_1; l_2 u_2]
end
# Homogeneous response (mutant fitness jointly inferred)
function CI_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff, MLL)
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
    u_1 = find_zero(LL_ratio_1, (m_UT, m_max(mc_max_UT,eff)))
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
    u_2 = find_zero(LL_ratio_2, (m_S, m_max(mc_max_S,eff)))
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
    u_3 = Inf
    try
        u_3 = find_zero(LL_ratio_3, (inv_fit_m, Inf))
    catch err
    end
    function min_M(P)
        if -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], P[2], P[3], eff) - MLL > chisq_1_95/2
            return Inf
        else
            return P[2]/P[1]
        end
    end
    res = Optim.optimize(min_M, [m_UT, m_S, inv_fit_m])
    l_4 = Optim.minimum(res)
    function max_M(P)
        if -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, P[1], P[2], P[3], eff) - MLL > chisq_1_95/2
            return Inf
        else
            return -P[2]/P[1]
        end
    end
    res = Optim.optimize(max_M, [m_UT, m_S, inv_fit_m])
    u_4 = -Optim.minimum(res)
    return [l_1 u_1; l_2 u_2; l_3 u_3; l_4 u_4]
end
# Heterogeneous response (fixed fraction and rel. division rate of on-cells)
function CI_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m, S, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, eff, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P[1], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
            res = Optim.optimize(log_likelihood_1, [S])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P1, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max_S,eff)))
    function LL_ratio_2(para)
        if para == S
            return -chisq_1_95/2
        else 
            log_likelihood_2(P) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], para, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
            res = Optim.optimize(log_likelihood_2, [m])
            P1 = Optim.minimizer(res)[1]
            return -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, para, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL - chisq_1_95/2
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(LL_ratio_2, (0., S))
    catch err
    end
    u_2 = find_zero(LL_ratio_2, (S, Inf))
    function min_m_on(P)
        if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
            return Inf
        else
            return P[2]*P[1]
        end
    end
    res = Optim.optimize(min_m_on, [m, S])
    l_3 = Optim.minimum(res)
    function max_m_on(P)
        if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
            return Inf
        else
            return -P[2]*P[1]
        end
    end
    res = Optim.optimize(max_m_on, [m, S])
    u_3 = -Optim.minimum(res)
    return [l_1 u_1; l_2 u_2; l_3 u_3]
end
# Heterogeneous response (fraction of on-cells fixed, rel. division rate of on-cells inferred)
function CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, inv_fit_m, eff, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P[1], f_on, P[2], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_1, [S, rel_div_on])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P1, f_on, P2, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max_S,eff)))
    function LL_ratio_2(para) 
        if para == S
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], para, f_on, P[2], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_2, [m, rel_div_on])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, para, f_on, P2, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(LL_ratio_2, (0., S))
    catch err
    end
    u_2 = find_zero(LL_ratio_2, (S, Inf))
    function LL_ratio_3(para)
        if para == rel_div_on
            return -chisq_1_95/2
        elseif para == 0.
            log_likelihood_r0(P) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
            res = Optim.optimize(log_likelihood_r0, [m, S])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL - chisq_1_95/2
        else 
            log_likelihood_3(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], f_on, para, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_3, [m, S])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, P2, f_on, para, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_3 = 0.
    try
        l_3 = find_zero(LL_ratio_3, (0., rel_div_on))
    catch err
    end
    u_3 = find_zero(LL_ratio_3, (rel_div_on, Inf))
    function min_m_on(P)
        if P[3] == 0. 
            if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                return Inf
            else 
                return P[2]*P[1]
            end
        else
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], f_on, P[3], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return P[2]*P[1]
            end
        end
    end
    res = Optim.optimize(min_m_on, [m, S, rel_div_on])
    l_4 = Optim.minimum(res)
    function max_m_on(P)
        if P[3] == 0. 
            if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                return Inf
            else 
                return -P[2]*P[1]
            end
        else
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], f_on, P[3], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return -P[2]*P[1]
            end
        end
    end
    res = Optim.optimize(max_m_on, [m, S, rel_div_on])
    u_4 = -Optim.minimum(res)
    return [l_1 u_1; l_2 u_2; l_3 u_3; l_4 u_4]
end
# Heterogeneous response (fraction of on-cells inferred, rel. division rate of on-cells fixed or inferred)
function CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, infer_r::Bool, inv_fit_m, eff, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P[1], P[2], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_1, [S, f_on])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, P1, P2, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = 0.
    try
        l_1 = find_zero(LL_ratio_1, (0., m))
    catch err
    end
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max_S,eff)))
    function LL_ratio_2(para) 
        if para == S
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], para, P[2], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_2, [m, f_on])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, para, P2, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_2 = 0.
    try 
        l_2 = find_zero(LL_ratio_2, (0., S))
    catch err
    end
    u_2 = find_zero(LL_ratio_2, (S, Inf))
    function LL_ratio_3(para) 
        if para == f_on
            return -chisq_1_95/2
        else
            log_likelihood_3(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], para, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
            res = Optim.optimize(log_likelihood_3, [m, S])
            P1, P2 = Optim.minimizer(res)
            return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, P2, para, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
        end
    end
    l_3 = 0.
    try 
        l_3 = find_zero(LL_ratio_3, (0., f_on))
    catch err
    end
    u_3 = 1.
    try 
        u_3 = find_zero(LL_ratio_3, (f_on, 1.))
    catch err
    end
    if infer_r == false
        function min_m_on(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return P[2]*P[1]*(1-P[3])/P[3]
            end
        end
        res = Optim.optimize(min_m_on, [m, S, f_on])
        l_4 = Optim.minimum(res)
        function max_m_on(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return -P[2]*P[1]*(1-P[3])/P[3]
            end
        end
        res = Optim.optimize(max_m_on, [m, S, f_on])
        u_4 = -Optim.minimum(res)
        function min_mu_inc(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return P[2]*(1-P[3])/P[3]
            end
        end
        res = Optim.optimize(min_mu_inc, [m, S, f_on])
        l_5 = Optim.minimum(res)
        function max_mu_inc(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return -P[2]*(1-P[3])/P[3]
            end
        end
        res = Optim.optimize(max_mu_inc, [m, S, f_on])
        u_5 = -Optim.minimum(res)
        function min_M(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return (1-P[3])*(1+P[2])
            end
        end
        res = Optim.optimize(min_M, [m, S, f_on])
        l_6 = Optim.minimum(res)
        function max_M(P)
            if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                return Inf
            else
                return -(1-P[3])*(1+P[2])
            end
        end
        res = Optim.optimize(max_M, [m, S, f_on])
        u_6 = -Optim.minimum(res)
        return [l_1 u_1; l_2 u_2; l_3 u_3; l_4 u_4; l_5 u_5; l_6 u_6]
    else
        function LL_ratio_r(para)
            if para == rel_div_on
                return -chisq_1_95/2
            elseif para == 0.
                log_likelihood_r0(P) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
                res = Optim.optimize(log_likelihood_r0, [m, S])
                P1, P2 = Optim.minimizer(res)
                return -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL - chisq_1_95/2
            else 
                log_likelihood_r(P) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], para, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
                res = Optim.optimize(log_likelihood_r, [m, S, f_on])
                P1, P2, P3 = Optim.minimizer(res)
                return -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P1, P2, P3, para, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL - chisq_1_95/2
            end
        end
        l_r = 0.
        try
            l_r = find_zero(LL_ratio_r, (0., rel_div_on))
        catch err
        end
        u_r = find_zero(LL_ratio_r, (rel_div_on, Inf))
        function min_r_m_on(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return P[2]*P[1]
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return P[2]*P[1]
                end
            end
        end
        res = Optim.optimize(min_r_m_on, [m, S, f_on, rel_div_on])
        l_4 = Optim.minimum(res)
        if l_4 == Inf
            l_4 = 0.
        end
        function max_r_m_on(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return -P[2]*P[1]
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return -P[2]*P[1]
                end
            end
        end
        res = Optim.optimize(max_r_m_on, [m, S, f_on, rel_div_on])
        u_4 = -Optim.minimum(res)
        if u_4 == -Inf
            u_4 = Inf
        end
        function min_r_mu_inc(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return P[2]*(1-P[3])/P[3]
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return P[2]*(1-P[3])/P[3]
                end
            end
        end
        res = Optim.optimize(min_r_mu_inc, [m, S, f_on, rel_div_on])
        l_5 = Optim.minimum(res)
        if l_5 == Inf
            l_5 = 0.
        end
        function max_r_mu_inc(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return -P[2]*(1-P[3])/P[3]
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return -P[2]*(1-P[3])/P[3]
                end
            end
        end
        res = Optim.optimize(max_r_mu_inc, [m, S, f_on, rel_div_on])
        u_5 = -Optim.minimum(res)
        if u_5 == -Inf
            u_5 = Inf
        end
        function min_r_M(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return (1-P[3])*(1+P[2])
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return (1-P[3])*(1+P[2])
                end
            end
        end
        res = Optim.optimize(min_r_M, [m, S, f_on, rel_div_on])
        l_6 = Optim.minimum(res)
        if l_6 == Inf
            l_6 = 0.
        end
        function max_r_M(P)
            if P[4] == 0. 
                if -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on) - MLL > chisq_1_95/2
                    return Inf
                else 
                    return -(1-P[3])*(1+P[2])
                end
            else
                if -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, P[1], P[2], P[3], P[4], q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff) - MLL > chisq_1_95/2
                    return Inf
                else
                    return -(1-P[3])*(1+P[2])
                end
            end
        end
        res = Optim.optimize(max_r_mu_inc, [m, S, f_on, rel_div_on])
        u_6 = -Optim.minimum(res)
        if u_6 == -Inf
            u_6 = Inf
        end
        return [l_1 u_1; l_2 u_2; l_3 u_3; l_r u_r; l_4 u_4; l_5 u_5; l_6 u_6]
    end
end

m_max(mc_max, eff::Bool) = 10*mc_max
m_max(mc_max, eff::Union{Float64,Tuple{Float64,Bool}}) = 10*mc_max/eff[1]
m_max(mc_max, eff) = 10*mc_max/min(eff[1][1], eff[2][1])