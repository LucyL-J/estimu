# Mutant fitness fixed -> q0, q coeffs. can be calculated outside of LL function
# Single fluctuation assay
function log_likelihood_m(mc_counts, mc_max, m, q0, q)
    if m <= 0.
        return -Inf
    else
        p = mudi_K(mc_max, m, q0, q)
        ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Mutant fitness inferred -> q0, q coeffs. are calculated inside LL function, and this calculation depends on plating efficiency
# Single fluctuation assay
function log_likelihood_m_fitm(mc_counts, mc_max, m, inv_fit_m, eff::Bool)
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        # Plating efficiency = 1 
        q0 = -1
        q = q_coeffs(mc_max, inv_fit_m)
        p = mudi_K(mc_max, m, q0, q)
        ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
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
        # Plating efficiency < 1
        q0 = q0_coeff(inv_fit_m, eff[1])
        q = q_coeffs(mc_max, inv_fit_m, eff)
        p = mudi_K(mc_max, m, q0, q)
        ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Mutant fitness = 0 -> Mutant count distribution = Poisson(m*eff)
function log_likelihood_m_fitm(mc_counts, mc_max, m_eff)
    if m_eff <= 0.
        return -Inf
    else
        p = mudi_K(mc_max, m_eff)
        ll = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Heterogeneous population with on-cells and off-cells
# Rel. division rate of on-cells = 0 -> q0, q coeffs. can be calculated outside of LL function
function log_likelihood_m_S(mc_counts, mc_max, m_off, S, q0_off, q_off, q0_on, q_on)
    if m_off <= 0. || S < 0.
        return -Inf
    else
        p = mudi_K(mc_max, m_off, q0_off, q_off, S*m_off, q0_on, q_on)
        ll  = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Rel. division rate of on-cells > 0 -> q0, q coeffs. for untreated calculated outside LL function, but for stressful inside and depending on plating efficiency
function log_likelihood_m_S_div_f(mc_counts, mc_max, m_off, S, f_on, rel_div_on, q0_off, q_off, inv_fit_m, eff::Bool)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        # Plating efficiency = 1
        q0_on = -1
        q_on = q_coeffs(mc_max, ifit)
        p = mudi_K(mc_max, m_off*scale_f(f_on, rel_div_on), q0_off, q_off, S*m_off*scale_f(f_on, rel_div_on), q0_on, q_on)
        ll  = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_S_div_f(mc_counts, mc_max, m_off, S, f_on, rel_div_on, q0_off, q_off, inv_fit_m, eff)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        # Plating efficiency < 1
        q0_on = q0_coeff(ifit, eff[1])
        q_on = q_coeffs(mc_max, ifit, eff)
        p = mudi_K(mc_max, m_off*scale_f(f_on, rel_div_on), q0_off, q_off, S*m_off*scale_f(f_on, rel_div_on), q0_on, q_on)
        ll  = sum(mc_counts .* log.(p))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end

# Pair of fluctuation assays: untreated and stressful
function log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m, q0_UT, q_UT, q0_S, q_S)
    if m <= 0.
        return -Inf
    else
        p_UT = mudi_K(mc_max_UT, m, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m*N_ratio, q0_S, q_S)
        ll = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end 
end
# Plating efficiency same for untreated and stressful -> q0, q coeffs. only calculated once and then subsampled
function log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, inv_fit_m, eff::Bool)
    # Number of mutations m, and inverse mutant fitness inv_fit_m jointly inferred each
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        # Plating efficiency = 1
        q0 = -1
        q = q_coeffs(mc_max, inv_fit_m)
        @views p_UT = mudi_K(mc_max_UT, m, q0, q[1:mc_max_UT])
        @views p_S = mudi_K(mc_max_S, m*N_ratio, q0, q[1:mc_max_S])
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff::Bool)
    # Number of mutations for untreated m_UT and stressful m_S inferred separately
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        # Plating efficiency = 1
        q0 = -1
        q = q_coeffs(mc_max, inv_fit_m)
        @views p_UT = mudi_K(mc_max_UT, m_UT, q0, q[1:mc_max_UT])
        @views p_S = mudi_K(mc_max_S, m_S, q0, q[1:mc_max_S])
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, inv_fit_m, eff::Union{Float64,Tuple{Float64,Bool}})
    # Number of mutations m jointly inferred
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        # Plating efficiency < 1
        q0 = q0_coeff(inv_fit_m, eff[1])
        q = q_coeffs(mc_max, inv_fit_m, eff)
        @views p_UT = mudi_K(mc_max_UT, m, q0, q[1:mc_max_UT])
        @views p_S = mudi_K(mc_max_S, m*N_ratio, q0, q[1:mc_max_S])
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff::Union{Float64,Tuple{Float64,Bool}})
    # Number of mutations for untreated m_UT and stressful m_S inferred separately
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        # Plating efficiency < 1
        q0 = q0_coeff(inv_fit_m, eff[1])
        q = q_coeffs(mc_max, inv_fit_m, eff)
        @views p_UT = mudi_K(mc_max_UT, m_UT, q0, q[1:mc_max_UT])
        @views p_S = mudi_K(mc_max_S, m_S, q0, q[1:mc_max_S])
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Mutant fitness untreated/stressful inferred separately -> q0, q coeffs. need to be calculated separately
function log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, x, N_ratio, m, inv_fit_m_UT, inv_fit_m_S, eff::Bool)
    if m <= 0. || inv_fit_m_UT <= 0. || inv_fit_m_S <= 0.
        return -Inf
    else
        # Plating efficiency = 1
        q0 = -1
        q_UT = q_coeffs(mc_max_UT, inv_fit_m_UT)
        q_S = q_coeffs(mc_max_S, inv_fit_m_S)
        p_UT = mudi_K(mc_max_UT, m, q0, q_UT)
        p_S = mudi_K(mc_max_S, m*N_ratio, q0, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, x, N_ratio, m, inv_fit_m_UT, inv_fit_m_S, eff::Union{Float64,Tuple{Float64,Bool}})
    if m <= 0. || inv_fit_m_UT <= 0. || inv_fit_m_S <= 0.
        return -Inf
    else
        # Plating efficiency < 1
        q0_UT = q0_coeff(inv_fit_m_UT, eff[1])
        q0_S = q0_coeff(inv_fit_m_S, eff[1])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m_UT, eff)
        q_S = q_coeffs(mc_max_S, inv_fit_m_S, eff)
        p_UT = mudi_K(mc_max_UT, m, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m*N_ratio, q0_S, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, x, m, inv_fit_m_UT, inv_fit_m_S, eff)
    if m <= 0. || inv_fit_m_UT <= 0. || inv_fit_m_S <= 0.
        return -Inf
    else
        # Plating efficiency different untreated/stressful
        q0_UT = q0_coeff(inv_fit_m_UT, eff[1][1])
        q0_S = q0_coeff(inv_fit_m_S, eff[2][1])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m_UT, eff[1])
        q_S = q_coeffs(mc_max_S, inv_fit_m_S, eff[2])
        p_UT = mudi_K(mc_max_UT, m, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m*N_ratio, q0_S, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Plating efficiency different for untreated and stressful -> q0, q coeffs. need to be calculated separately
function log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, x, m, inv_fit_m, eff)
    if m <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1][1])
        q0_S = q0_coeff(inv_fit_m, eff[2][1])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1])
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[2])
        p_UT = mudi_K(mc_max_UT, m, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m*N_ratio, q0_S, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, x, m_UT, m_S, inv_fit_m, eff)
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1][1])
        q0_S = q0_coeff(inv_fit_m, eff[2][1])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1])
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[2])
        p_UT = mudi_K(mc_max_UT, m_UT, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m_S, q0_S, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Mutant fitnes = 0 -> Mutant count distributions = Poisson(m_UT*eff_UT) and Poisson(m_S*eff_S)
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_eff_UT, m_eff_S)
    if m_eff_UT <= 0. || m_eff_S <= 0.
        return -Inf
    else
        p_UT = mudi_K(mc_max_UT, m_eff_UT)
        p_S = mudi_K(mc_max_S, m_eff_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end

# Heterogeneous population with on-cells and off-cells
# Rel. division rate of on-cells = 0 -> q0, q coeffs. can be calculated outside of LL function
function log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
    if m_off <= 0. || S < 0.
        return -Inf
    else
        p_UT = mudi_K(mc_max_UT, m_off, q0_UT, q_UT)
        p_S = mudi_K(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
# Rel. division rate of on-cells > 0 -> q0, q coeffs. for untreated calculated outside LL function, but for stressful inside and depending on plating efficiency
function log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff::Bool)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        p_UT = mudi_K(mc_max_UT, m_off, q0_UT, q_UT)
        N_ratio *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        # Plating efficiency = 1
        q0_S_on = -1
        q_S_on = q_coeffs(mc_max_S, ifit)
        p_S = mudi_K(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        p_UT = mudi_K(mc_max_UT, m_off, q0_UT, q_UT)
        N_ratio *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        # Plating efficiency < 1
        q0_S_on = q0_coeff(ifit, eff[1])
        q_S_on = q_coeffs(mc_max_S, ifit, eff)
        p_S = mudi_K(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 0.
            return ll
        else
            return -Inf
        end
    end
end

# Probability generating function method used to set initial value of the mutation-supply ratio S for the maximum likelihood estimation, based on
# Gillet-Markowska, A., Louvel, G., & Fischer, G. (2015). bz-rates: A web tool to estimate mutation rates from fluctuation analysis. G3: Genes, Genomes, Genetics, 5(11), 2323–2327. https://doi.org/10.1534/g3.115.019836
function empirical_pgf(z, x) # Empirical probability generating function calculated from observed data
    g = 0
    for i in x
        g += z^i
    end
    g /= length(x)
    return g
end
function initial_S(z, mc, m)             # Estimate the mutation-supply ratio for given z   
    if z == 0.
        return -log(empirical_pgf(z, mc))/m - 1
    else
        return -log(empirical_pgf(z, mc))/(m*(1-z)) + log(1-z)/z
    end
end
function initial_S(mc, m, z_values::Int) # Estimate the mutation-supply ratio by averaging over a number of z values 
    S = 0.
    for i = 0:z_values-1
        S += initial_S(i/z_values, mc, m)
    end
    if S == Inf || S < 0.
        return 0.
    else
        return S/z_values
    end
end
function initial_f(mc, N_ratio, Nf_S, m, S, rel_div_on)
    if m == 0.
        return 0.
    else
        mu_inc = max(1.,median(mc))/(m*N_ratio)
        f_upper = 1 - mu_inc/(S+1)                                                          
        if f_upper <= 0.
            f_upper = 1/mu_inc
        end
        f_lower = - log(1-f_upper) / log(Nf_S)                                                                                                                                                                                                                                        
        f_on = f_lower/(1 - rel_div_on)                                                                         
        return minimum([maximum([f_lower, f_on]), f_upper])
    end
end