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
    q0, q = coeffs(mc_max, 1/fit_m, eff)
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

# Mutation rate estimation from pair of fluctuation assays under permissive/stressful cond. without change in mutation rate
# Input
# mc_UT: Mutant counts untreated
# Nf_UT: Average final population size untreated
# mc_S: Mutant counts under stressful cond.
# Nf_S: Average final population size under stressful cond.
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 but can be set to different value(s) if known from separate experiment(s)
#        If only one value is given (instead of a vector), mutant fitness is constrained to be equal under permissive/stressful cond.
# cond_S: Condition, by default = "S" for stressful 

# Mutant fitness fixed in the inference
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{Float64}, fit_m::Vector{Float64}=[1., 1.]; cond_S="S") 
    if fit_m[1] ==  fit_m[2] == 1.                                         
        m = "No SIM"
    else                              
        m = "No SIM (diff. mutant fitness)"
    end
    cond_m = "UT+"*cond_S
    N_ratio = Nf_S/Nf_UT
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutant fitness"], condition=[cond_m, "UT", cond_S], status=["jointly inferred", "set to input", "set to input"])       
    msel_res = DataFrame(model=[m], status=["-"])   
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    mc_max = max(mc_max_UT, mc_max_S)
    if eff[1] == eff[2] && fit_m[1] == fit_m[2]
        q0, q = coeffs(mc_max, 1/fit_m[1], eff[1])
        q0_UT = q0
        q0_S = q0
        q_UT = q[1:mc_max_UT]
        q_S = q[1:mc_max_S]
    else
        q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
        q0_S, q_S = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    end
    # 1 inference parameter: Number of mutations under untreated+stressful cond.                                            
    LL(para) = -log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, q0_UT, q_UT, q0_S, q_S)
    # Maximum mutant count observed overall, used as an upper bound for the inference parameter
    res = Optim.optimize(LL, 0., mc_max)                                      
    if Optim.converged(res) == true
        est_res.MLE = [Optim.minimizer(res)/Nf_UT, fit_m[1], fit_m[2]]           
        msel_res.LL = [-Optim.minimum(res)]
        msel_res.AIC = [2 + 2*Optim.minimum(res)]
        # Number of data points = total number of parallel cultures
        msel_res.BIC = [1*log(length(mc_UT)+length(mc_S)) + 2*Optim.minimum(res)] 
    else
        est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end
    return est_res, msel_res
end