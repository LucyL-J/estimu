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
# mc_S:  Mutant counts under stressful cond.
# Nf_S:  Average final population size under stressful cond.
# eff:   Plating efficiency 
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 for untreated and stressful cond. but can be set to different value(s) if known from separate experiment(s)
# cond_S: Condition, by default = "S" for stressful 

# Mutant fitness fixed in the inference
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Vector{Float64}=[1., 1.]; cond_S="S") 
    if fit_m[1] ==  fit_m[2] == 1.                                         
        m = "No SIM"
    else                              
        m = "No SIM (diff. mutant fitness)"
    end
    N_ratio = Nf_S/Nf_UT
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutant fitness"], condition=["UT+"*cond_S, "UT", cond_S], status=["jointly inferred", "set to input", "set to input"])       
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

# Mutation rate estimation from pair of fluctuation assays under permissive/stressful cond. using the homogeneous-response model
# Input
# mc_UT: Mutant counts untreated
# Nf_UT: Average final population size untreated
# mc_S:  Mutant counts under stressful cond.
# Nf_S:  Average final population size under stressful cond.
# eff:   Plating efficiency 
# Optional
# fit_m: Mutant fitness
#        By default fit_m=1 for untreated and stressful cond. but can be set to different value(s) if known from separate experiment(s)
#        Unless set as a joint inference parameter via fitm=false, inference under permissive/stressful cond(s). is independent
# cond_S: Condition, by default = "S" for stressful 

# Mutant fitness under permissive/stressful cond(s). independent (either fixed in the inference, or inferred separately)
function estimu_hom(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m=[1., 1.]; cond_s="S")
    if typeof(fit_m) == Vector{Float64} 
        if fit_m[1] == fit_m[2] == 1.                                         
            m = "Homogeneous"
        elseif fit_m[1] == fit_m[2]
            m = "Homogeneous (constr. fitness)"
        else
            m = "Homogeneous (unconstr. mutant fitness)"
        end
    else
        m = "Homogeneous (unconstr. mutant fitness)"
    end
    # Estimation for permissive cond.
	est_res_p, msel_res_p = estimu(mc_UT, Nf_UT, eff[1], fit_m[1])
    if msel_res_p.LL[1] != -Inf
        # Estimation for stressful cond.
        est_res_s, msel_res_s = estimu(mc_S, Nf_S, eff[2], fit_m[2], cond=cond_s)
        if msel_res_s.LL[1] != -Inf
            est_res_p = vcat(est_res_p, est_res_s)
            if typeof(fit_m[2]) ==  Bool
                s = "calc. from 2&4"
            else
                s = "set to input"
            end
            push!(est_res_p, ["Ratio mutant fitness", cond_s*"/UT", s, est_res_s.MLE[2]/est_res_p.MLE[2]])
            push!(est_res_p, ["Fold change mutation rate", cond_s*"/UT", "calc. from 1&3", est_res_s.MLE[1]/est_res_p.MLE[1]])
            msel_res_p.LL += msel_res_s.LL
            msel_res_p.AIC += msel_res_s.AIC
        else
            push!(est_res_p, ["Mutation rate", cond_s, "failed", 0.])
            push!(est_res_p, ["Mutant fitness", cond_s, "failed", -1.])
            msel_res_p.LL = [-Inf]
            msel_res_p.AIC = [Inf]
        end
    end
    msel_res_p.BIC = [sum((typeof(fit_m[1])==Bool)+(sum(typeof(fit_m[2])==Bool))) * log(length(mc_UT)+length(mc_S)) - 2*msel_res_p.LL[1]]
    msel_res_p.model = [m]      
    return est_res_p, msel_res_p
end

# Mutant fitness jointly inferred (constrained to be equal under permissive/stressful cond(s).)
function estimu_hom(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Bool; cond_s="S")
    parameter=["Mutation rate", "Mutant fitness", "Mutation rate", "Mutant fitness", "Ratio mutant fitness", "Fold change mutation rate"]
    condition = ["UT", "UT+"*cond_s, cond_s, "UT+"*cond_s, "", cond_s*"/UT"]
    status = ["inferred", "jointly inferred", "inferred", "jointly inferred", "constr.", "calc. from 1&3"]
    est_res = DataFrame(parameter=parameter, condition=condition, status=status)
    msel_res = DataFrame(model=["Homogeneous (constr. mutant fitness)"], status=["-"])                              
    # 3 inference parameters: Number of mutations under permissive/stressful cond., mutant fitness
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    mc_max = max(mc_max_UT, mc_max_S)
    if eff[1] == eff[2] == 1
        LL_eff_1(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3])
        res = Optim.optimize(LL_eff_1, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
    elseif eff[1] == eff[2]
        if eff[1] < 0.5
            LL_small_eff(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3], eff[1], true)
            res = Optim.optimize(LL_small_eff, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
        else
            LL(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3], eff[1])
            res = Optim.optimize(LL, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
        end
    else
        if eff[1] < 0.5
            if eff[2] < 0.5
                LL_small_eff_12(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, para[1], para[2], para[3], eff, true)
                res = Optim.optimize(LL_small_eff_12, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
            else
                LL_small_eff_1(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, para[1], para[2], para[3], (eff[1], true, eff[2]))
                res = Optim.optimize(LL_small_eff_1, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
            end
        else
            if eff[2] < 0.5
                LL_small_eff_2(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, para[1], para[2], para[3], (eff[1], eff[2], true))
                res = Optim.optimize(LL_small_eff_2, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
            else
                LL_12(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3], eff)
                res = Optim.optimize(LL_12, [initial_m(mc_UT, 1000), initial_m(mc_S, 1000), 1.]) 
            end
        end     
    end
	if Optim.converged(res) == true
		p = Optim.minimizer(res)
        est_res.MLE = [p[1]/Nf_UT, p[3], p[2]/Nf_S, p[3], 1., p[2]/p[1] * Nf_UT/Nf_S]
        msel_res.LL = [-Optim.minimum(res)]
        msel_res.AIC = [2*(length(Nf_S)+2) + 2*Optim.minimum(res)]
        msel_res.BIC = [log(length(mc_UT)+sum([length(mc) for mc in mc_S]))*(length(Nf_S)+2) + 2*Optim.minimum(res)] 
	else
		est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
	end
	return est_res, msel_res
end
