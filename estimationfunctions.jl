include("mutantcountdistributions.jl")
include("loglikelihoods.jl")
using StatsBase, DataFrames, Optim

# Mutation rate estimation algorithms
# Return two data frames: 
#       (i) Maximum likelihood estimates and 95% confidence intervals 
#       (ii) Log-likelihood and AIC+BIC values

# Mutation rate estimation from single fluctuation assay using the standard model with optional differential mutant fitness
# Input
# mc:  Mutant counts
# Nf:  Average final population size
# eff: Plating efficiency
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 but can be set to different value if known from separate experiment
#        Alternatively, mutant fitness is inferred if fit_m=false
# cond: Condition, by default = "UT" for untreated
function estimu(mc::Vector{Int}, Nf, eff, fit_m::Float64=1.; cond="UT")
    # Create the output data frames
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness"], condition=[cond, cond])
    msel_res = DataFrame()
    # Model depends on the value of the mutant fitness
    if fit_m == 1.
        msel_res = DataFrame(model=["Standard"])
    else
        msel_res = DataFrame(model=["Standard (diff. mutant fitness)"])
    end       
    msel_res.status = ["-"]       
    # Pre-inference calculations
    mc_max = maximum(mc)
    mc_counts = counts(mc, 0:mc_max)
    if eff == 1.
        q0 = -1
        if fit_m == 1.
            q = q_coeffs(mc_max)
        else
            q = q_coeffs(mc_max, 1/fit_m)
        end
    else
        q0 = q0_coeff(1/fit_m, eff)
        if eff < 0.5
            q = q_coeffs(mc_max, 1/fit_m, eff, true)
        else
            q = q_coeffs(mc_max, 1/fit_m, eff)
        end
    end
    # 1 inference parameter: Number of mutations 
    LL(para) = -log_likelihood_m(mc_counts, mc_max, para, q0, q)
    res = Optim.optimize(LL, 0., mc_max)
    if Optim.converged(res) == true
        est_res.status = ["inferred", "set to input"]
        # Mutation rate is calculated from m and the final population size
        est_res.MLE = [Optim.minimizer(res)/Nf, fit_m]     
        msel_res.LL = [-Optim.minimum(res)]
        msel_res.AIC = [2 + 2*Optim.minimum(res)]
        msel_res.BIC = [log(length(mc)) + 2*Optim.minimum(res)]                         
    else
        est_res.status = fill("failed", length(est_res.parameter))                                                      
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end 
    return est_res, msel_res
end
# Mutant fitness not given -> inferred 
function estimu(mc::Vector{Int}, Nf, eff, fit_m::Bool; cond="UT")
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness"], condition=[cond, cond])
    msel_res = DataFrame(model=["Standard (diff. mutant fitness)"], status=["-"])                                                        
	mc_max = maximum(mc)
    mc_counts = counts(mc, 0:mc_max)
    # 2 inference parameters: Number of mutations, mutant fitness
    if eff == 1
        LL_eff_1(para) = -log_likelihood_m_fitm(mc_counts, mc_max, para[1], para[2])
        res = Optim.optimize(LL_eff_1, [initial_m(mc, 1000), 1.]) 
    elseif eff < 0.5
        LL_small_eff(para) = -log_likelihood_m_fitm(mc_counts, mc_max, para[1], para[2], eff, true)
        res = Optim.optimize(LL_small_eff, [initial_m(mc, 1000), 1.]) 
    else
        LL(para) = -log_likelihood_m_fitm(mc_counts, mc_max, para[1], para[2], eff)
        res = Optim.optimize(LL, [initial_m(mc, 1000), 1.]) 
    end
	if Optim.converged(res) == true
		p = Optim.minimizer(res)
		est_res.status = ["inferred", "inferred"]
		est_res.MLE = [p[1]/Nf, 1/p[2]]                         
        msel_res.LL = [-Optim.minimum(res)]
        msel_res.AIC = [4 + 2*Optim.minimum(res)]         
        msel_res.BIC = [2*log(length(mc)) + 2*Optim.minimum(res)]                       
	else
        est_res.status = fill("failed", length(est_res.parameter))
		msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
	end
	return est_res, msel_res                                                
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
function initial_mu_sup(z, mc, m)             # Estimate the mutation-supply ratio for given z   
    if z == 0.
        return -log(empirical_pgf(z, mc))/m - 1
    else
        return -log(empirical_pgf(z, mc))/(m*(1-z)) + log(1-z)/z
    end
end
function initial_mu_sup(mc, m, z_values::Int) # Estimate the mutation-supply ratio by averaging over a number of z values 
    mu_sup = 0.
    for i = 0:z_values-1
        mu_sup += initial_mu_sup(i/z_values, mc, m)
    end
    if mu_sup == Inf || mu_sup < 0.
        return 0.
    else
        return mu_sup/z_values
    end
end