include("mutantcountdistributions.jl")
include("loglikelihoods.jl")
include("confidenceintervals.jl")
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
    # Model depends on the value of the mutant fitness
    if fit_m == 1.
        M = "Standard"
    else
        M = "Standard (diff. mutant fitness)"
    end     
    # Create the output data frames
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness"], condition=[cond, cond])
    msel_res = DataFrame(model=[M], status=["-"])   
    # Pre-inference calculations
    mc_max = maximum(mc)
    mc_counts = counts(mc, 0:mc_max)
    q0, q = coeffs(mc_max, 1/fit_m, eff)
    # 1 inference parameter: Number of mutations 
    LL(para) = -log_likelihood_m(mc_counts, mc_max, para, q0, q)
    res = Optim.optimize(LL, 0., mc_max, iterations=10^4)
    if Optim.converged(res) == true
        est_res.status = ["inferred", "set to input"]
        m = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        # Mutation rate is calculated from m and the final population size
        est_res.MLE = [m/Nf, fit_m]     
        try
            b = CI_m(mc_counts, mc_max, m, q0, q, eff, MLL)
            est_res.lower_bound = [b[1]/Nf, fit_m]
            est_res.upper_bound = [b[2]/Nf, fit_m] 
        catch
            est_res.lower_bound = [0., fit_m]
            est_res.upper_bound = [Inf, fit_m] 
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [2 + 2*MLL]
        msel_res.BIC = [log(length(mc)) + 2*MLL]                         
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
    # Different cases regarding partial plating
    if eff == 1
        eff = false
    elseif eff < 0.5
        eff = (eff, true)
    end
    # 2 inference parameters: Number of mutations, mutant fitness
    LL(para) = -log_likelihood_m_fitm(mc_counts, mc_max, para[1], para[2], eff)
    res = Optim.optimize(LL, [max(1.,median(mc)), 1.], iterations=10^4) 
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.status = ["inferred", "inferred"]
		est_res.MLE = [p[1]/Nf, 1/p[2]]   
        try
            b = CI_m_fitm(mc_counts, mc_max, p[1], p[2], eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf, 1/b[2,2]]
            est_res.upper_bound = [b[1,2]/Nf, 1/b[2,1]]
        catch
            est_res.lower_bound = [0., 0.]
            est_res.upper_bound = [Inf, Inf]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [4 + 2*MLL]         
        msel_res.BIC = [2*log(length(mc)) + 2*MLL]                       
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
        M = "No SIM"
    else                              
        M = "No SIM (diff. mutant fitness)"
    end
    N_ratio = Nf_S/Nf_UT
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutant fitness"], condition=["UT+"*cond_S, "UT", cond_S])       
    msel_res = DataFrame(model=[M], status=["-"])   
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
    res = Optim.optimize(LL, 0., mc_max, iterations=10^4)                                      
    if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "set to input", "set to input"]
        m = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [m/Nf_UT, fit_m[1], fit_m[2]]   
        try
            b = CI_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, q0_UT, q_UT, q0_S, q_S, eff, MLL)        
            est_res.lower_bound = [b[1]/Nf_UT, fit_m[1], fit_m[2]]
            est_res.upper_bound = [b[2]/Nf_UT, fit_m[1], fit_m[2]]
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2]]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2]]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [2 + 2*MLL]
        # Number of data points = total number of parallel cultures
        msel_res.BIC = [1*log(length(mc_UT)+length(mc_S)) + 2*MLL] 
    else
        est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end
    return est_res, msel_res
end
# Mutant fitness jointly inferred
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Bool; cond_S="S") 
    N_ratio = Nf_S/Nf_UT
    est_res = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutant fitness"], condition=["UT+"*cond_S, "UT+"*cond_S, "UT+"*cond_S])       
    msel_res = DataFrame(model=["No SIM (diff. mutant fitness)"], status=["-"])                            
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    mc_max = max(mc_max_UT, mc_max_S)
    if eff[1] == eff[2] == 1
        eff = false
    elseif eff[1] == eff[2]
        if eff[1] < 0.5
            eff = (eff[1], true)
        else
            eff = eff[1]
        end
    else
        eff_UT, eff_S = eff
        if eff[1] < 0.5
            eff_UT = (eff[1], true)
        end
        if eff[2] < 0.5
            eff_S = (eff[2], true)
        end
        eff = (eff_UT, eff_S)
    end
    # 2 inference parameters: Number of mutations, mutant fitness
    LL(para) = -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para[1], para[2], eff)
    res = Optim.optimize(LL, [max(1.,median(mc_UT),median(mc_S)), 1.], iterations=10^4) 
	if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "jointly inferred", "jointly inferred"]
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, 1/p[2], 1/p[2]]   
        try
            b = CI_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, p[1], p[2], eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, 1/b[2,2], 1/b[2,2]]
            est_res.upper_bound = [b[1,2]/Nf_UT, 1/b[2,1], 1/b[2,1]]
        catch
            est_res.lower_bound = [0., 0., 0.]
            est_res.upper_bound = [Inf, Inf, Inf]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [4 + 2*MLL]
        msel_res.BIC = [2*log(length(mc_UT)+length(mc_S)) + 2*MLL] 
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
function estimu_hom(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Union{Vector{Float64},Tuple{Bool,Bool},BitVector}=[1., 1.]; cond_S="S")
    if typeof(fit_m) == Vector{Float64} 
        if fit_m[1] == fit_m[2] == 1.                                         
            M = "Homogeneous"
        elseif fit_m[1] == fit_m[2]
            M = "Homogeneous (constr. fitness)"
        else
            M = "Homogeneous (unconstr. mutant fitness)"
        end
    else
        M = "Homogeneous (unconstr. mutant fitness)"
    end
    # Estimation for permissive cond.
	est_res_UT, msel_res_UT = estimu(mc_UT, Nf_UT, eff[1], fit_m[1])
    if msel_res_UT.LL[1] != -Inf
        # Estimation for stressful cond.
        est_res_S, msel_res_S = estimu(mc_S, Nf_S, eff[2], fit_m[2], cond=cond_S)
        if msel_res_S.LL[1] != -Inf
            mc_max_UT = maximum(mc_UT)
            mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
            mc_max_S = maximum(mc_S)
            mc_counts_S = counts(mc_S, 0:mc_max_S)
            if typeof(fit_m) == Vector{Float64}
                q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
                q0_S, q_S = coeffs(mc_max_S, 1/fit_m[2], eff[2])
                b_M = CI_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, est_res_UT.MLE[1]*Nf_UT, est_res_S.MLE[1]*Nf_S, q0_UT, q_UT, q0_S, q_S, -msel_res_UT.LL[1]-msel_res_S.LL[1])
                b = [b_M; fit_m[1]/fit_m[2] fit_m[1]/fit_m[2]]
            else
                b = CI_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, est_res_UT.MLE[1]*Nf_UT, est_res_S.MLE[1]*Nf_S, 1/est_res_UT.MLE[2], 1/est_res_S.MLE[2], eff, -msel_res_UT.LL[1]-msel_res_S.LL[1])
            end
            est_res_UT = vcat(est_res_UT, est_res_S)
            if typeof(fit_m[2]) ==  Bool
                s = "calc. from 2&4"
            else
                s = "set to input"
            end
            push!(est_res_UT, ["Ratio mutant fitness", cond_S*"/UT", s, est_res_S.MLE[2]/est_res_UT.MLE[2], 1/b[2,2], 1/b[2,1]])
            push!(est_res_UT, ["Fold change mutation rate", cond_S*"/UT", "calc. from 1&3", est_res_S.MLE[1]/est_res_UT.MLE[1], b[1,1]*Nf_UT/Nf_S, b[1,2]*Nf_UT/Nf_S])
            msel_res_UT.LL += msel_res_S.LL
            msel_res_UT.AIC += msel_res_S.AIC
        else
            push!(est_res_UT, ["Mutation rate", cond_S, "failed", 0., 0., 0.])
            push!(est_res_UT, ["Mutant fitness", cond_S, "failed", -1., -1., -1.])
            msel_res_UT.LL = [-Inf]
            msel_res_UT.AIC = [Inf]
        end
    end
    msel_res_UT.BIC = [sum((typeof(fit_m[1])==Bool)+(sum(typeof(fit_m[2])==Bool))) * log(length(mc_UT)+length(mc_S)) - 2*msel_res_UT.LL[1]]
    msel_res_UT.model = [M]      
    return est_res_UT, msel_res_UT
end
# Mutant fitness jointly inferred (constrained to be equal under permissive/stressful cond(s).)
function estimu_hom(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Bool; cond_S="S")
    parameter=["Mutation rate", "Mutant fitness", "Mutation rate", "Mutant fitness", "Ratio mutant fitness", "Fold change mutation rate"]
    condition = ["UT", "UT+"*cond_S, cond_S, "UT+"*cond_S, "", cond_S*"/UT"]
    status = ["inferred", "jointly inferred", "inferred", "jointly inferred", "constr.", "calc. from 1&3"]
    est_res = DataFrame(parameter=parameter, condition=condition, status=status)
    msel_res = DataFrame(model=["Homogeneous (constr. mutant fitness)"], status=["-"])                              
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    mc_max = max(mc_max_UT, mc_max_S)
    if eff[1] == eff[2] == 1
        eff = false
    elseif eff[1] == eff[2]
        if eff[1] < 0.5
            eff = (eff[1], true)
        else
            eff = eff[1]
        end
    else
        eff_UT, eff_S = eff
        if eff[1] < 0.5
            eff_UT = (eff[1], true)
        end
        if eff[2] < 0.5
            eff_S = (eff[2], true)
        end
        eff = (eff_UT, eff_S)
    end
    # 3 inference parameters: Number of mutations under permissive/stressful cond., mutant fitness
    LL(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3], eff)
    res = Optim.optimize(LL, [max(1.,median(mc_UT)), max(1.,median(mc_S)), 1.], iterations=10^4) 
	if Optim.converged(res) == true
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, 1/p[3], p[2]/Nf_S, 1/p[3], 1., p[2]/p[1] * Nf_UT/Nf_S]
        try
            b = CI_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, p[1], p[2], p[3], eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, 1/b[3,2], b[2,1]/Nf_S, 1/b[3,2], 1., b[4,1] * Nf_UT/Nf_S]
            est_res.upper_bound = [b[1,2]/Nf_UT, 1/b[3,1], b[2,2]/Nf_S, 1/b[3,1], 1., b[4,2] * Nf_UT/Nf_S]
        catch
            est_res.lower_bound = [0., 0., 0., 0., 1., 0.]
            est_res.upper_bound = [Inf, Inf, Inf, Inf, 1., Inf]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [6 + 2*MLL]
        msel_res.BIC = [3*log(length(mc_UT)+length(mc_S)) + 2*MLL] 
	else
		est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
	end
	return est_res, msel_res
end

# Input
# mc_UT: Mutant counts untreated
# Nf_UT: Average final population size untreated
# mc_S:  Mutant counts under stressful cond.
# Nf_S:  Average final population size under stressful cond.
# eff:   Plating efficiency 
# Optional
# f_on: Fraction of on-cells
#       Either inferred or set to different value(s) if known from a separate experiment
# rel_div_on: Relative division rate of on-cells compared to off-cells
#             By default set to rel_div_on=0., inferred if rel_div_on=false
#             For rel_div_on=0., the fraction of on-cells cannot be inferred
# cond_S: Condition, by default = "S" for stressful 

# Fraction and relative division rate of on-cells given and fixed in the inference
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Float64, rel_div_on::Float64=0., fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    est_res = DataFrame(parameter=["Mutation rate off-cells", "Mutant fitness", "Mutant fitness", "Mutation-supply ratio", "Mutation rate on-cells", "Fraction on-cells", "Rel. division rate on-cells", "Rel. mutation rate on-cells", "Fold change mean mutation rate"])
	est_res.condition = [["UT+"*cond_S, "UT"]; fill(cond_S, 6); cond_S*"/UT"]
    if rel_div_on == 0.
        M = ["Heterogeneous (zero division rate on-cells)"]
	else
        M = ["Heterogeneous"]
	end
    msel_res = DataFrame(model=M, status=["-"])
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    N_ratio = Nf_S/Nf_UT
    # Calculate initial values
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    if rel_div_on == 0.
        q0_S_on = -eff[2]
        q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    else
        sf = scale_f(f_on, rel_div_on)
        N_ratio *= sf
        ifit = inverse_fit_on(f_on, rel_div_on)/fit_m[2]
        q0_S_on, q_S_on = coeffs(mc_max_S, ifit, eff[2])
    end
    # 2 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio
    LL(para) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
    res = Optim.optimize(LL, [m, S], iterations=10^4)                                          
    if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "set to input", "set to input", "inferred", "calc. from 1,4&6", "set to input", "set to input", "calc. from 4&6", "calc. from 4&6"] 
        p = Optim.minimizer(res)     
        MLL = Optim.minimum(res)                                                          
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[2]*p[1]*(1-f_on)/(f_on*Nf_UT), f_on, rel_div_on, p[2]*(1-f_on)/f_on, (1-f_on)*(1+p[2])]
        try
            b = CI_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[3,1]*(1-f_on)/(f_on*Nf_UT), f_on, rel_div_on, b[2,1]*(1-f_on)/f_on, (1-f_on)*(1+b[2,1])]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[3,2]*(1-f_on)/(f_on*Nf_UT), f_on, rel_div_on, b[2,2]*(1-f_on)/f_on, (1-f_on)*(1+b[2,2])]   
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., f_on, rel_div_on, 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, f_on, rel_div_on, Inf, Inf]  
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [4 + 2*MLL]         
        msel_res.BIC = [2*log(length(mc_UT)+length(mc_S)) + 2*MLL]   
    else
        est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end   
    return est_res, msel_res
end
# Fraction of on-cells fixed, rel. division rate of on-cells inferred
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Float64, rel_div_on::Bool, fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    est_res = DataFrame(parameter=["Mutation rate off-cells", "Mutant fitness", "Mutant fitness", "Mutation-supply ratio", "Mutation rate on-cells", "Fraction on-cells", "Rel. division rate on-cells", "Rel. mutation rate on-cells", "Fold change mean mutation rate"])
	est_res.condition = [["UT+"*cond_S, "UT"]; fill(cond_S, 6); cond_S*"/UT"]
    msel_res = DataFrame(model=["Heterogeneous"])
    msel_res.status = ["-"] 
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    N_ratio = Nf_S/Nf_UT
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2]) 
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    if eff[2] == 1
        eff = false
    elseif eff[2] < 0.5
        eff = (eff[2], true)
    else
        eff = eff[2]
    end
    # 3 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, rel. division rate on-cells                              
    LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], f_on, para[3], q0_UT, q_UT, q0_S_off, q_S_off, 1/fit_m[2], eff)
    res = Optim.optimize(LL, [m, S, 1.], iterations=10^4)                                    
    if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "set to input", "set to input", "inferred", "calc. from 1,4&6", "set to input", "inferred", "calc. from 4&6", "calc. from 4&6"] 
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[2]*p[1]*(1-f_on)/(f_on*Nf_UT), f_on, p[3], p[2]*(1-f_on)/f_on, (1-f_on)*(1+p[2])]
        try
            b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], f_on, p[3], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, 1/fit_m[2], eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[4,1]*(1-f_on)/(f_on*Nf_UT), f_on, b[3,1], b[2,1]*(1-f_on)/f_on, (1-f_on)*(1+b[2,1])]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[4,2]*(1-f_on)/(f_on*Nf_UT), f_on, b[3,2], b[2,2]*(1-f_on)/f_on, (1-f_on)*(1+b[2,2])] 
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., f_on, 0., 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, f_on, Inf, Inf, Inf]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [6 + 2*MLL]         
        msel_res.BIC = [3*log(length(mc_UT)+length(mc_S)) + 2*MLL] 
    else
        est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end
    return est_res, msel_res
end
# Fraction of on-cells not given; inferred if the rel. division rate of on-cells is non-zero
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Bool, rel_div_on::Float64=0., fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    parameter = ["Mutation rate off-cells", "Mutant fitness", "Mutant fitness"]
    condition = ["UT+"*cond_S, "UT", cond_S]
    status = ["jointly inferred", "set to input", "set to input"]
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    N_ratio = Nf_S/Nf_UT
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    # For zero rel. division rate on-cells -> Fraction of on-cells cannot be inferred
    if rel_div_on == 0.
        parameter = [parameter; ["Mutation-supply ratio", "Rel. division rate on-cells"]]
        condition = [condition; [cond_S, cond_S]]
        status = [status; ["inferred", "set to input"]]
        est_res = DataFrame(parameter=parameter, condition=condition, status=status)
        msel_res = DataFrame(model=["Heterogeneous (zero division rate on-cells)"], status=["-"])
        # 2 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio  
        LL_0(para) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
        res = Optim.optimize(LL_0, [m, S], iterations=10^4)                     
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            MLL = Optim.minimum(res)
            est_res.MLE  = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], 0.]
            try
                b = CI_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, eff, MLL)
                est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], 0.]
                est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], 0.] 
            catch
                est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0.]
                est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, 0.]

            end
            msel_res.LL = [-MLL]
            msel_res.AIC = [4 + 2*MLL]         
            msel_res.BIC = [2*log(length(mc_UT)+length(mc_S)) + 2*MLL]
        else
            est_res.status = fill("failed", length(est_res.parameter))
            msel_res.LL = [-Inf]
            msel_res.AIC = [Inf]
            msel_res.BIC = [Inf]
        end
    # For non-zero rel. division rate on-cells -> Fraction of on-cells inferred
    else
        parameter = [parameter; ["Mutation-supply ratio", "Mutation rate on-cells", "Fraction on-cells", "Rel. division rate on-cells", "Rel. mutation rate on-cells", "Fold change mean mutation rate"]]
        condition = [condition; [fill(cond_S, 5); cond_S*"/UT"]]
        status = [status; ["inferred", "calc. from 1,4&6", "inferred", "set to input", "calc. from 4&6", "calc. from 4&6"]]
        # Calculate the initial value for optimisation
        f_on = initial_f(mc_S, N_ratio, Nf_S, m, S, rel_div_on)
        est_res = DataFrame(parameter=parameter, condition=condition, status=status)
        msel_res = DataFrame(model=["Heterogeneous"], status=["-"])
        if eff[2] == 1
            eff = false
        elseif eff[2] < 0.5
            eff = (eff[2], true)
        else
            eff = eff[2]
        end
        # 3 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, fraction of on-cells                              
        LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], para[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, 1/fit_m[2], eff)
        res = Optim.optimize(LL, [m, S, f_on], iterations=10^4) 
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            MLL = Optim.minimum(res)
            est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[2]*p[1]*(1-p[3])/(p[3]*Nf_UT), p[3], rel_div_on, p[2]*(1-p[3])/p[3], (1-p[3])*(1+p[2])]
            try
                b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], p[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, false, 1/fit_m[2], eff, MLL)
                est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[4,1]/Nf_UT, b[3,1], rel_div_on, b[5,1], b[6,1]]
                est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[4,2]/Nf_UT, b[3,2], rel_div_on, b[5,2], b[6,2]] 
            catch
                est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., 0., rel_div_on, 0., 0.]
                est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, 1., rel_div_on, Inf, Inf]

            end
            msel_res.LL = [-MLL]
            msel_res.AIC = [6 + 2*MLL]         
            msel_res.BIC = [3*log(length(mc_UT)+length(mc_S)) + 2*MLL] 
        else
            est_res.status = fill("failed", length(est_res.parameter))
            msel_res.LL = [-Inf]
            msel_res.AIC = [Inf]
            msel_res.BIC = [Inf]
        end   
    end
    return est_res, msel_res 
end
# Relative division rate on-cells and fraction of on-cells not given -> both inferred
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Bool, rel_div_on::Bool, fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    est_res = DataFrame(parameter=["Mutation rate off-cells", "Mutant fitness", "Mutant fitness", "Mutation-supply ratio", "Mutation rate on-cells", "Fraction on-cells", "Rel. division rate on-cells", "Rel. mutation rate on-cells", "Fold change mean mutation rate"])
	est_res.condition = [["UT+"*cond_S, "UT"]; fill(cond_S, 6); cond_S*"/UT"]
    msel_res = DataFrame(model=["Heterogeneous"], status=["-"])
    mc_max_UT = maximum(mc_UT)
    mc_counts_UT = counts(mc_UT, 0:mc_max_UT)
    mc_max_S = maximum(mc_S)
    mc_counts_S = counts(mc_S, 0:mc_max_S)
    N_ratio = Nf_S/Nf_UT
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    f_on = initial_f(mc_S, N_ratio, Nf_S, m, S, 0.)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    if eff[2] == 1
        eff = false
    elseif eff[2] < 0.5
        eff = (eff[2], true)
    else
        eff = eff[2]
    end
    # 4 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, relative division rate on-cells, fraction of on-cells
    LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], para[3], para[4], q0_UT, q_UT, q0_S_off, q_S_off, fit_m[2], eff)
    res = Optim.optimize(LL, [m, S, f_on, 1.], iterations=10^4)                                                         
    if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "set to input", "set to input", "inferred", "calc. from 1,4&6", "inferred", "inferred", "calc. from 4&6", "calc. from 4&6"] 
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[2]*p[1]*(1-p[3])/(p[3]*Nf_UT), p[3], p[4], p[2]*(1-p[3])/p[3], (1-p[3])*(1+p[2])]                                                          
        try
            b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], p[3], p[4], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, true, 1/fit_m[2], eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[5,1]/Nf_UT, b[3,1], b[4,1], b[6,1], b[7,1]]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[5,2]/Nf_UT, b[3,2], b[4,2], b[6,2], b[7,2]] 
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., 0., 0., 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, 1., Inf, Inf, Inf]
        end
        msel_res.LL = [-MLL]
        msel_res.AIC = [8 + 2*MLL]         
        msel_res.BIC = [4*log(length(mc_UT)+length(mc_S)) + 2*MLL]  
    else
        est_res.status = fill("failed", length(est_res.parameter))
        msel_res.LL = [-Inf]
        msel_res.AIC = [Inf]
        msel_res.BIC = [Inf]
    end
    return est_res, msel_res
end