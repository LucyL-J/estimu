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
# Without change in mutation rate (mutant fitness inferred separately)
function CI_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, inv_fit_m_UT, inv_fit_m_S, eff, MLL)
    function LL_ratio_1(para)
        if para == m
            return -chisq_1_95/2
        else
            log_likelihood_1(P) = -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para, P[1], P[2], eff)
            res = Optim.optimize(log_likelihood_1, [inv_fit_m_UT, inv_fit_m_S])
            P_res = Optim.minimizer(res)
            return -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para, P_res[1], P_res[2], eff) - MLL - chisq_1_95/2
        end
    end
    l_1 = find_zero(LL_ratio_1, (0., m))
    u_1 = find_zero(LL_ratio_1, (m, m_max(mc_max,eff)))
    function LL_ratio_2(para)
        if para == inv_fit_m_UT
            return -chisq_1_95/2
        else
            log_likelihood_2(P) = -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P[1], para, P[2], eff)
            res = Optim.optimize(log_likelihood_2, [m, inv_fit_m_S])
            P_res = Optim.minimizer(res)
            return -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P_res[1], para, P_res[2], eff) - MLL - chisq_1_95/2
        end
    end
    l_2 = find_zero(LL_ratio_2, (0., inv_fit_m_UT))
    u_2 = Inf
    try
        u_2 = find_zero(LL_ratio_2, (inv_fit_m_UT, Inf))
    catch err
    end
    function LL_ratio_3(para)
        if para == inv_fit_m_S
            return -chisq_1_95/2
        else
            log_likelihood_3(P) = -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P[1], P[2], para, eff)
            res = Optim.optimize(log_likelihood_3, [m, inv_fit_m_UT])
            P_res = Optim.minimizer(res)
            return -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, P_res[1], P_res[2], para, eff) - MLL - chisq_1_95/2
        end
    end
    l_3 = find_zero(LL_ratio_3, (0., inv_fit_m_S))
    u_3 = Inf
    try
        u_3 = find_zero(LL_ratio_3, (inv_fit_m_S, Inf))
    catch err
    end
    l_4, u_4 = CI_fitm(inv_fit_m_UT, inv_fit_m_S, l_2, l_3, u_2, u_3)
    return [l_1 u_1; l_2 u_2; l_3 u_3; l_4 u_4]
end
CI_fitm(inv_fit_m_UT, inv_fit_m_S, min_inv_fit_m_UT, min_inv_fit_m_S, max_inv_fit_m_UT, max_inv_fit_m_S) = [min(inv_fit_m_S/max_inv_fit_m_UT, min_inv_fit_m_S/inv_fit_m_UT), max(max_inv_fit_m_S/inv_fit_m_UT, inv_fit_m_S/min_inv_fit_m_UT)]
# Homogeneous response (fixed mutant fitness)
CI_m(m_UT, m_S, min_m_UT, min_m_S, max_m_UT, max_m_S) = [min(m_S/max_m_UT, min_m_S/m_UT) max(max_m_S/m_UT, m_S/min_m_UT)]
# Homogeneous response (mutant fitness inferred separately)
CI_m_fitm(m_UT, m_S, min_m_UT, min_m_S, max_m_UT, max_m_S, inv_fit_m_UT, inv_fit_m_S, min_inv_fit_m_UT, min_inv_fit_m_S, max_inv_fit_m_UT, max_inv_fit_m_S) = [min(m_S/max_m_UT, min_m_S/m_UT) max(max_m_S/m_UT, m_S/min_m_UT); min(inv_fit_m_S/max_inv_fit_m_UT, min_inv_fit_m_S/inv_fit_m_UT) max(max_inv_fit_m_S/inv_fit_m_UT, inv_fit_m_S/min_inv_fit_m_UT)]
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
    return [l_1 u_1; l_2 u_2; l_3 u_3; min(l_2/m_UT, m_S/u_1) max(u_2/m_UT, m_S/l_1)]
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
    return [l_1 u_1; l_2 u_2; min(l_2*m, l_1*S) max(u_2*m, u_1*S)]
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
    return [l_1 u_1; l_2 u_2; l_3 u_3; min(l_2*m, l_1*S) max(u_2*m, u_1*S)]
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
        return [l_1 u_1; l_2 u_2; l_3 u_3; min(l_2*m*(1-f_on)/f_on, S*l_1*(1-f_on)/f_on, S*m*(1-u_3)/u_3) max(u_2*m*(1-f_on)/f_on, S*u_1*(1-f_on)/f_on, S*m*(1-l_3)/l_3); min(l_2*(1-f_on)/f_on, S*(1-u_3)/u_3) max(u_2*(1-f_on)/f_on, S*(1-l_3)/l_3); min((1-f_on)*(1+l_2), (1-u_3)*(1+S)) max((1-f_on)*(1*u_2), (1-l_3)*(1+S))]
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
        return [l_1 u_1; l_2 u_2; l_r u_r; l_3 u_3; min(l_2*m*(1-f_on)/f_on, S*l_1*(1-f_on)/f_on, S*m*(1-u_3)/u_3) max(u_2*m*(1-f_on)/f_on, S*u_1*(1-f_on)/f_on, S*m*(1-l_3)/l_3); min(l_2*(1-f_on)/f_on, S*(1-u_3)/u_3) max(u_2*(1-f_on)/f_on, S*(1-l_3)/l_3); min((1-f_on)*(1+l_2), (1-u_3)*(1+S)) max((1-f_on)*(1*u_2), (1-l_3)*(1+S))]
    end
end

m_max(mc_max, eff::Bool) = 10*mc_max
m_max(mc_max, eff::Union{Float64,Tuple{Float64,Bool}}) = 10*mc_max/eff[1]
m_max(mc_max, eff) = 10*mc_max/min(eff[1][1], eff[2][1])