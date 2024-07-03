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
function log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m, q0_UT, q_UT, q0_S, q_S)
    if m <= 0.
        return -Inf
    else
        p_UT = mudi(mc_max_UT, m, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m*N_ratio, q0_S, q_S)
        @views ll = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
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
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m)
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = -1
        q = q_coeffs(mc_max, inv_fit_m)
        @views p_UT = mudi(mc_max_UT, m_UT, q0, q[1:mc_max_UT])
        @views p_S = mudi(mc_max_S, m_S, q0, q[1:mc_max_S])
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff::Float64)
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = q0_coeff(inv_fit_m, eff)
        q = q_coeffs(mc_max, inv_fit_m, eff)
        @views p_UT = mudi(mc_max_UT, m_UT, q0, q[1:mc_max_UT])
        @views p_S = mudi(mc_max_S, m_S, q0, q[1:mc_max_S])
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, m_UT, m_S, inv_fit_m, eff::Float64, small_eff::Bool)
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0 = q0_coeff(inv_fit_m, eff)
        q = q_coeffs(mc_max, inv_fit_m, eff, true)
        @views p_UT = mudi(mc_max_UT, m_UT, q0, q[1:mc_max_UT])
        @views p_S = mudi(mc_max_S, m_S, q0, q[1:mc_max_S])
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, inv_fit_m, eff::Vector{Float64}, small_eff::Bool)
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1])
        q0_S = q0_coeff(inv_fit_m, eff[2])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1], true)
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[2], true)
        p_UT = mudi(mc_max_UT, m_UT, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m_S, q0_S, q_S)
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, inv_fit_m, eff::Tuple{Float64,Bool,Float64})
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1])
        q0_S = q0_coeff(inv_fit_m, eff[3])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1], true)
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[3])
        p_UT = mudi(mc_max_UT, m_UT, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m_S, q0_S, q_S)
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, inv_fit_m, eff::Tuple{Float64,Float64,Bool})
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1])
        q0_S = q0_coeff(inv_fit_m, eff[2])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1])
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[2], true)
        p_UT = mudi(mc_max_UT, m_UT, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m_S, q0_S, q_S)
        @views ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, m_UT, m_S, inv_fit_m, eff::Vector{Float64})
    if m_UT <= 0. || m_S <= 0. || inv_fit_m <= 0.
        return -Inf
    else
        q0_UT = q0_coeff(inv_fit_m, eff[1])
        q0_S = q0_coeff(inv_fit_m, eff[2])
        q_UT = q_coeffs(mc_max_UT, inv_fit_m, eff[1])
        q_S = q_coeffs(mc_max_S, inv_fit_m, eff[2])
        p_UT = mudi(mc_max_UT, m_UT, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m_S, q0_S, q_S)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end

function log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
    if m_off <= 0. || S < 0.
        return -Inf
    else
        p_UT = mudi(mc_max_UT, m_off, q0_UT, q_UT)
        p_S = mudi(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        p_UT = mudi(mc_max_UT, m_off, q0_UT, q_UT)
        N_ratio *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        q0_S_on = -1
        q_S_on = q_coeffs(mc_max_S, ifit)
        p_S = mudi(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
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
        p_UT = mudi(mc_max_UT, m_off, q0_UT, q_UT)
        N_ratio *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        q0_S_on = q0_coeff(ifit, eff)
        q_S_on = q_coeffs(mc_max_S, ifit, eff)
        p_S = mudi(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end
function log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, m_off, S, f_on, rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, inv_fit_m, eff, small_eff::Bool)
    if m_off <= 0. || S < 0. || rel_div_on <= 0. || f_on < 0. || f_on >= 1. 
        return -Inf
    else
        p_UT = mudi(mc_max_UT, m_off, q0_UT, q_UT)
        N_ratio *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)*inv_fit_m
        q0_S_on = q0_coeff(ifit, eff)
        q_S_on = q_coeffs(mc_max_S, ifit, eff, true)
        p_S = mudi(mc_max_S, m_off*N_ratio, q0_S_off, q_S_off, S*m_off*N_ratio, q0_S_on, q_S_on)
        ll  = sum(mc_counts_UT .* log.(p_UT)) + sum(mc_counts_S .* log.(p_S))
        if !isnan(ll) && ll < 1.
            return ll
        else
            return -Inf
        end
    end
end

# Probability generating function method used to set initial values of the maximum likelihood estimation, based on
# Gillet-Markowska, A., Louvel, G., & Fischer, G. (2015). bz-rates: A web tool to estimate mutation rates from fluctuation analysis. G3: Genes, Genomes, Genetics, 5(11), 2323–2327. https://doi.org/10.1534/g3.115.019836
function empirical_pgf(z, x) # Empirical probability generating function calculated from observed data
    g = 0
    for i in x
        g += z^i
    end
    g /= length(x)
    return g
end
function initial_m(z, mc)                     # Estimate number of mutations for given z        
    if z == 0.
        return log(empirical_pgf(z, mc)) 
    else
        return z/((1-z)*log(1-z)) * log(empirical_pgf(z, mc))
    end
end
function initial_m(mc, z_values::Int)         # Estimate number of mutations by averaging over a number of z values
    m = 0.
    for i = 0:z_values-1
        m += initial_m(i/z_values, mc)
    end
    return max(m/z_values, 0.)
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
function initial_f(mc, N_ratio, Nf_S, m, S, rel_div_on, z_values::Int)
    mu_inc = initial_m(mc, z_values)/(m*N_ratio)
    f_upper = 1 - mu_inc/(S+1)                                                          
    if f_upper <= 0.
        f_upper = 1/mu_inc
    end
    f_lower = - log(1-f_upper) / log(Nf_S)                                                                                                                                                                                                                                        
    f_on = f_lower/(1 - rel_div_on)                                                                         
    return minimum([maximum([f_lower, f_on]), f_upper])
end