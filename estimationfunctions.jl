include("mutantcountdistributions.jl")
include("loglikelihoods.jl")
include("confidenceintervals.jl")
include("goodnessoffit.jl")
using StatsBase, DataFrames, Optim

# Mutation rate estimation algorithms
# Return two data frames: 
#       (i) Maximum likelihood estimates and 95% confidence intervals 
#       (ii) Log-likelihood and AIC, AIC_corrected + BIC values

# Global estimation meta parameters
R_gof = 10^4            # Number of replicates for the goodness-of-fit test
max_mc_cutoff = 1000    # Maximum mutant count at which the mutant count distribution is cutoff for the goodness-of-fit test

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

# To create the output data frames with default/failed values
est_res_standard(cond) = DataFrame(parameter=["Mutation rate", "Mutant fitness"], condition=[cond, cond], status=["inferred", "set to input"], MLE=[0.,1.], lower_bound=[0.,1.], upper_bound=[Inf,1.])
msel_res_standard(M, cond) = DataFrame(model=[M], condition=[cond], AIC=[Inf], AIC_corr=[Inf], BIC=[Inf], LL=[-Inf], p_value=[0.], cutoff=[0], tail_prob=[1.], calc_time=[0.]) 

function estimu(mc::Vector{Int}, Nf, eff, fit_m::Float64=1.; cond="UT")
    start_time = time()
    est_res = est_res_standard(cond)
    # Model depends on the value of the mutant fitness
    M = (fit_m == 1.) ? "Standard" : "Standard (diff. mutant fitness)"     
    msel_res = msel_res_standard(M, cond)
    # Pre-inference calculations
    mc_max, mc_counts, num_c = extract_mc(mc)
    q0, q = coeffs(mc_max, 1/fit_m, eff)
    # 1 inference parameter: Number of mutations 
    num_para = 1
    LL(para) = -log_likelihood_m(mc_counts, mc_max, para, q0, q)
    res = Optim.optimize(LL, 0., mc_max, iterations=10^4)
    if Optim.converged(res) == true
        m = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        # Mutation rate is calculated from m and the final population size
        est_res.MLE = [m/Nf, fit_m]     
        try # Confidence interval calculation
            b = CI_m(mc_counts, mc_max, m, q0, q, eff, MLL)
            est_res.lower_bound = [b[1]/Nf, fit_m]
            est_res.upper_bound = [b[2]/Nf, fit_m] 
        catch # If the confidence interval calculation fails, set the bounds to 0 and Inf
            est_res.lower_bound = [0., fit_m]
            est_res.upper_bound = [Inf, fit_m] 
        end
        # Calculate the selection criteria: AIC, AIC corrected (for small sample size), BIC
        msel_res[1,3:5] = selection_crit(num_para, num_c, MLL)    
        msel_res.LL = [-MLL]
        # Goodness-of-fit test: distribution of log-likelihoods of randomly drawn mutant counts under the estimated parameters
        LLs, mc_cutoff, p_cutoff = LL_dist(num_c, Nf, m/Nf, fit_m, eff) # Also returns the cutoff value and tail probability in the pmf calculation
        msel_res.p_value = [1 - ecdf(LLs)(MLL)]
        msel_res.cutoff = [mc_cutoff]
        msel_res.tail_prob = [p_cutoff] 
        end_time = time()
        msel_res.calc_time = [end_time - start_time]  
    else
        est_res.status .= "failed"
    end 
    return est_res, msel_res, LLs
end
# Mutant fitness not given -> inferred 
function estimu(mc::Vector{Int}, Nf, eff, fit_m::Bool; cond="UT")
    start_time = time()
    est_res = est_res_standard(cond)
    msel_res = msel_res_standard("Standard (diff. mutant fitness)", cond)
	mc_max, mc_counts, num_c = extract_mc(mc)
    # Different cases regarding partial plating
    eff_conv = eff
    if eff == 1
        eff_conv = false
    elseif eff < 0.5
        eff_conv = (eff, true)
    end
    # 2 inference parameters: Number of mutations, mutant fitness
    num_para = 2
    LL(para) = -log_likelihood_m_fitm(mc_counts, mc_max, para[1], para[2], eff_conv)
    res = Optim.optimize(LL, [max(1.,median(mc)), 1.], iterations=10^4) 
    if Optim.converged(res) == true
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.status[2] = "inferred"
		est_res.MLE = [p[1]/Nf, 1/p[2]]   
        try
            b = CI_m_fitm(mc_counts, mc_max, p[1], p[2], eff_conv, MLL)
            est_res.lower_bound = [b[1,1]/Nf, 1/b[2,2]]
            est_res.upper_bound = [b[1,2]/Nf, 1/b[2,1]]
        catch
            est_res.lower_bound = [0., 0.]
            est_res.upper_bound = [Inf, Inf]
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c, MLL)   
        msel_res.LL = [-MLL]
        LLs, mc_cutoff, p_cutoff = LL_dist(num_c, Nf, p[1]/Nf, 1/p[2], eff)
        msel_res.p_value = [1 - ecdf(LLs)(MLL)]  
        msel_res.cutoff = [mc_cutoff]
        msel_res.tail_prob = [p_cutoff]
        end_time = time()
        msel_res.calc_time = [end_time - start_time]               
	else
        est_res.status .= "failed"
        LLs = zeros(Float64, R_gof) # Return empty LLs if the inference fails; needed for estimu_hom with separate estimation under permissive/stressful cond(s)
	end
	return est_res, msel_res, LLs                                                
end 

# Mutation rate estimation from pair of fluctuation assays under permissive/stressful cond.
msel_res_joint(M, cond_S) = DataFrame(model=[M,M,M], condition=["UT+"*cond_S,"UT",cond_S], AIC=[Inf,NaN,NaN], AIC_corr=[Inf,NaN,NaN], BIC=[Inf,NaN,NaN], LL=[-Inf,-Inf,-Inf], p_value=[0.,0.,0.], cutoff=[-1,0,0], tail_prob=[NaN,1.,1.], calc_time=[0.,NaN,NaN])   

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
est_res_joint_0(cond_S) = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutant fitness", "Ratio mutant fitness"], condition=["UT+"*cond_S, "UT", cond_S, cond_S*"/UT"], 
    status=["jointly inferred", "set to input", "set to input", "set to input"], MLE=[0., 1., 1., 1.], lower_bound=[0., 1., 1., 1.], upper_bound=[Inf, 1., 1., 1.]) 

# Mutant fitness fixed in the inference
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Vector{Float64}=[1., 1.]; cond_S="S") 
    start_time = time()
    est_res = est_res_joint_0(cond_S)
    M = (fit_m[1] == fit_m[2] == 1.) ? "No SIM" : (fit_m[1] == fit_m[2]) ? "No SIM (constr. mutant fitness)" : "No SIM (unconstr. mutant fitness)"
    msel_res = msel_res_joint(M, cond_S)
    N_ratio = Nf_S/Nf_UT
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    mc_max = max(mc_max_UT, mc_max_S)
    if eff[1] == eff[2] && fit_m[1] == fit_m[2]
        q0, q = coeffs(mc_max, 1/fit_m[1], eff[1])
        q0_UT = q0_S = q0
        q_UT = q[1:mc_max_UT]
        q_S = q[1:mc_max_S]
    else
        q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
        q0_S, q_S = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    end
    # 1 inference parameter: Number of mutations under untreated+stressful cond.                                            
    num_para = 1
    LL(para) = -log_likelihood_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para, q0_UT, q_UT, q0_S, q_S)
    # Maximum mutant count observed overall, used as an upper bound for the inference parameter
    res = Optim.optimize(LL, 0., mc_max, iterations=10^4)                                      
    if Optim.converged(res) == true
        m = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [m/Nf_UT, fit_m[1], fit_m[2], fit_m[2]/fit_m[1]]   
        try
            b = CI_joint_m(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, m, q0_UT, q_UT, q0_S, q_S, eff, MLL)        
            est_res.lower_bound = [b[1]/Nf_UT, fit_m[1], fit_m[2], fit_m[2]/fit_m[1]]
            est_res.upper_bound = [b[2]/Nf_UT, fit_m[1], fit_m[2], fit_m[2]/fit_m[1]]
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], fit_m[2]/fit_m[1]]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], fit_m[2]/fit_m[1]]
        end
        # Sample size = total number of parallel cultures
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m(mc_counts_UT, mc_max_UT, m, q0_UT, q_UT)
        LL_S = log_likelihood_m(mc_counts_S, mc_max_S, m*N_ratio, q0_S, q_S)
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, m/Nf_UT, fit_m[1], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, m/Nf_UT, fit_m[2], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
    else
        est_res.status .= "failed"
    end
    return est_res, msel_res
end
# Mutant fitness jointly inferred
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Bool; cond_S="S") 
    start_time = time()
    est_res = est_res_joint_0(cond_S)
    msel_res = msel_res_joint("No SIM (constr. mutant fitness)", cond_S)
    N_ratio = Nf_S/Nf_UT  
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    mc_max = max(mc_max_UT, mc_max_S)
    eff_conv = convert_eff(eff)
    # 2 inference parameters: Number of mutations, mutant fitness
    num_para = 2
    LL(para) = -log_likelihood_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para[1], para[2], eff_conv)
    res = Optim.optimize(LL, [max(1.,median(mc_UT),median(mc_S)), 1.], iterations=10^4) 
	if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "jointly inferred", "jointly inferred", "constr."]
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, 1/p[2], 1/p[2], 1.]   
        try
            b = CI_joint_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, p[1], p[2], eff_conv, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, 1/b[2,2], 1/b[2,2], 1.]
            est_res.upper_bound = [b[1,2]/Nf_UT, 1/b[2,1], 1/b[2,1], 1.]
        catch
            est_res.lower_bound = [0., 0., 0., 1.]
            est_res.upper_bound = [Inf, Inf, Inf, 1.]
        end
        if eff[1] == eff[2]
            eff_conv = (eff_conv, eff_conv)
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, p[1], p[2], eff_conv[1])
        LL_S = log_likelihood_m_fitm(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], eff_conv[2])
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, 1/p[2], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, 1/p[2], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
	else
		est_res.status .= "failed"
	end
	return est_res, msel_res
end
# Mutant fitness inferred separately
function estimu_0(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Union{Tuple{Bool,Bool},BitVector}; cond_S="S") 
    start_time = time()
    est_res = est_res_joint_0(cond_S) 
    msel_res = msel_res_joint("No SIM (unconstr. mutant fitness)", cond_S)  
    N_ratio = Nf_S/Nf_UT                     
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    mc_max = max(mc_max_UT, mc_max_S)
    eff_conv = convert_eff(eff)
    # 3 inference parameters: Number of mutations, mutant fitness under untreated/stressful cond.
    num_para = 3
    LL(para) = -log_likelihood_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, para[1], para[2], para[3], eff_conv)
    res = Optim.optimize(LL, [max(1.,median(mc_UT),median(mc_S)), 1., 1.], iterations=10^4) 
	if Optim.converged(res) == true
        est_res.status = ["jointly inferred", "inferred", "inferred", "calculated from 2&3"]
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, 1/p[2], 1/p[3], p[2]/p[3]]
        CI_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, p[1], p[2], p[3], eff_conv, MLL)   
        try
            b = CI_joint_m_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, N_ratio, p[1], p[2], p[3], eff_conv, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, 1/b[2,2], 1/b[3,2], 1/b[4,2]]
            est_res.upper_bound = [b[1,2]/Nf_UT, 1/b[2,1], 1/b[3,1], 1/b[4,1]]
        catch
            est_res.lower_bound = [0., 0., 0., 0.]
            est_res.upper_bound = [Inf, Inf, Inf, Inf]
        end
        if eff[1] == eff[2]
            eff_conv = (eff_conv, eff_conv)
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, p[1], p[2], eff_conv[1])
        LL_S = log_likelihood_m_fitm(mc_counts_S, mc_max_S, p[1]*N_ratio, p[3], eff_conv[2])
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, 1/p[2], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, 1/p[3], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
	else
		est_res.status .= "failed"
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
    M = (typeof(fit_m) == Vector{Float64}) ? ((fit_m[1] == fit_m[2] == 1.) ? "Homogeneous" : ((fit_m[1] == fit_m[2]) ? "Homogeneous (constr. fitness)" : "Homogeneous (unconstr. mutant fitness)")) : "Homogeneous (unconstr. mutant fitness)"
    msel_res = msel_res_joint(M, cond_S)
    # Estimation for permissive cond.
	est_res_UT, msel_res_UT, LLs_UT = estimu(mc_UT, Nf_UT, eff[1], fit_m[1])
    msel_res[2,:] = msel_res_UT[1,:] # Copy the UT cond. to the second row of msel_res
    # Estimation for stressful cond.
    est_res_S, msel_res_S, LLs_S = estimu(mc_S, Nf_S, eff[2], fit_m[2], cond=cond_S)
    msel_res[3,:] = msel_res_S[1,:] # Copy the S cond. to the second row of msel_res
    est_res = vcat(est_res_UT, est_res_S)
    if msel_res_UT.LL[1] != -Inf && msel_res_S.LL[1] != -Inf
        start_time = time()
        mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
        if typeof(fit_m) == Vector{Float64}
            q0_S, q_S = coeffs(mc_max_S, 1/fit_m[2], eff[2])
            LL_S_UT = log_likelihood_m(mc_counts_S, mc_max_S, est_res_UT.MLE[1]*Nf_UT, q0_S, q_S)
            LLs_S_UT, mc_cutoff_S_UT, p_cutoff_S_UT = LL_dist(num_c_S, Nf_S, est_res_UT.MLE[1], fit_m[2], eff[2])
            b_M = CI_m(est_res_UT.MLE[1]*Nf_UT, est_res_S.MLE[1]*Nf_S, est_res_UT.lower_bound[1]*Nf_UT, est_res_S.lower_bound[1]*Nf_S, est_res_UT.upper_bound[1]*Nf_UT, est_res_S.upper_bound[1]*Nf_S)
            b = [b_M; fit_m[1]/fit_m[2] fit_m[1]/fit_m[2]]
        else
            LL_S_UT = log_likelihood_m_fitm(mc_counts_S, mc_max_S, est_res_UT.MLE[1]*Nf_UT, 1/est_res_UT.MLE[2], eff[2]) 
            LLs_S_UT, mc_cutoff_S_UT, p_cutoff_S_UT = LL_dist(num_c_S, Nf_S, est_res_UT.MLE[1], est_res_UT.MLE[2], eff[2])
            b = CI_m_fitm(est_res_UT.MLE[1]*Nf_UT, est_res_S.MLE[1]*Nf_S, est_res_UT.lower_bound[1]*Nf_UT, est_res_S.lower_bound[1]*Nf_S, est_res_UT.upper_bound[1]*Nf_UT, est_res_S.upper_bound[1]*Nf_S, 1/est_res_UT.MLE[2], 1/est_res_S.MLE[2], 1/est_res_UT.upper_bound[2], 1/est_res_S.upper_bound[2], 1/est_res_UT.lower_bound[2], 1/est_res_S.lower_bound[2])
        end
        s = (typeof(fit_m[2]) ==  Bool) ? "calc. from 2&4" : "set to input"
        push!(est_res, ["Ratio mutant fitness", cond_S*"/UT", s, est_res_S.MLE[2]/est_res_UT.MLE[2], 1/b[2,2], 1/b[2,1]])
        push!(est_res, ["Fold change mutation rate", cond_S*"/UT", "calc. from 1&3", est_res_S.MLE[1]/est_res_UT.MLE[1], b[1,1]*Nf_UT/Nf_S, b[1,2]*Nf_UT/Nf_S])
        num_para = 2 + sum((typeof(fit_m[1]) == Bool)+(sum(typeof(fit_m[2]) == Bool))) # Total number of inference parameters
        msel_res[1,3:5] = selection_crit(num_para, length(mc_UT)+length(mc_S), msel_res_UT.LL[1]+msel_res_S.LL[1])
        end_time = time()
        total_time = msel_res_UT.calc_time[1]+msel_res_S.calc_time[1] + (end_time-start_time)
        msel_res[1,6:end] = [msel_res_UT.LL[1]+msel_res_S.LL[1], 1 - ecdf(LLs_UT.+LLs_S)(-msel_res_UT.LL[1]-msel_res_S.LL[1]), -1, NaN, total_time]
        push!(msel_res, [M, "UT->"*cond_S, NaN, NaN, NaN, LL_S_UT, 1 - ecdf(LLs_S_UT)(-LL_S_UT), mc_cutoff_S_UT, p_cutoff_S_UT, total_time])
    else
        push!(msel_res, [M, "UT->"*cond_S, NaN, NaN, NaN, -Inf, 0., 0, 1., 0.])
    end
    msel_res.model .= M     
    return est_res, msel_res
end
# Mutant fitness jointly inferred (constrained to be equal under permissive/stressful cond(s).)
est_res_joint_hom(cond_S) = DataFrame(parameter=["Mutation rate", "Mutant fitness", "Mutation rate", "Mutant fitness", "Ratio mutant fitness", "Fold change mutation rate"], condition=["UT", "UT+"*cond_S, cond_S, "UT+"*cond_S, cond_S*"/UT", cond_S*"/UT"],
    status=["inferred", "jointly inferred", "inferred", "jointly inferred", "constr.", "calc. from 1&3"], MLE=[0., 1., 0., 1., 1., 1.], lower_bound=[0., 0., 0., 0., 1., 0.], upper_bound=[Inf, Inf, Inf, Inf, 1., Inf])

function estimu_hom(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, fit_m::Bool; cond_S="S")
    start_time = time()
    est_res = est_res_joint_hom(cond_S)
    msel_res = msel_res_joint("Homogeneous (constr. mutant fitness)", cond_S)
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    mc_max = max(mc_max_UT, mc_max_S)
    eff_conv = convert_eff(eff)
    # 3 inference parameters: Number of mutations under permissive/stressful cond., mutant fitness
    num_para = 3
    LL(para) = -log_likelihood_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, para[1], para[2], para[3], eff_conv)
    res = Optim.optimize(LL, [max(1.,median(mc_UT)), max(1.,median(mc_S)), 1.], iterations=10^4) 
	if Optim.converged(res) == true
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, 1/p[3], p[2]/Nf_S, 1/p[3], 1., p[2]/p[1] * Nf_UT/Nf_S]
        try
            b = CI_m_joint_fitm(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, mc_max, p[1], p[2], p[3], eff_conv, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, 1/b[3,2], b[2,1]/Nf_S, 1/b[3,2], 1., b[4,1] * Nf_UT/Nf_S]
            est_res.upper_bound = [b[1,2]/Nf_UT, 1/b[3,1], b[2,2]/Nf_S, 1/b[3,1], 1., b[4,2] * Nf_UT/Nf_S]
        catch
            est_res.lower_bound = [0., 0., 0., 0., 1., 0.]
            est_res.upper_bound = [Inf, Inf, Inf, Inf, 1., Inf]
        end
        if eff[1] == eff[2]
            eff_conv = (eff_conv, eff_conv)
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, p[1], p[3], eff_conv[1])
        LL_S = log_likelihood_m_fitm(mc_counts_S, mc_max_S, p[2], p[3], eff_conv[2])
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, 1/p[3], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[2]/Nf_S, 1/p[3], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
	else
		est_res.status .= "failed"
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
est_res_joint_het(cond_S) = DataFrame(parameter=["Mutation rate off-cells", "Mutant fitness", "Mutant fitness", "Mutation-supply ratio", "Rel. division rate on-cells", "Fraction on-cells", "Mutation rate on-cells", "Rel. mutation rate on-cells", "Fold change mean mutation rate"],
    condition=[["UT+"*cond_S, "UT"]; fill(cond_S, 6); cond_S*"/UT"], status=["jointly inferred", "set to input", "set to input", "inferred", "set to input", "set to input", "calc. from 1,4&6", "calc. from 4&6", "calc. from 4&6"],
    MLE=[0., 1., 1., 0., 0., 0., 0., 0., 0.], lower_bound=[0., 1., 1., 0., 0., 0., 0., 0., 0.], upper_bound=[Inf, 1., 1., Inf, Inf, Inf, Inf, Inf, Inf])

# Fraction and relative division rate of on-cells given and fixed in the inference
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Float64, rel_div_on::Float64=0., fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    start_time = time()
    est_res = est_res_joint_het(cond_S)
    M = (rel_div_on == 0.) ? "Heterogeneous (zero division rate on-cells)" : "Heterogeneous"
    msel_res = msel_res_joint(M, cond_S)
    N_ratio = Nf_S/Nf_UT
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
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
    num_para = 2
    LL(para) = -log_likelihood_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on)
    res = Optim.optimize(LL, [m, S], iterations=10^4)                                          
    if Optim.converged(res) == true 
        p = Optim.minimizer(res)     
        MLL = Optim.minimum(res)                                                          
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], rel_div_on, f_on, p[2]*p[1]*(1-f_on)/(f_on*Nf_UT), p[2]*(1-f_on)/f_on, (1-f_on)*(1+p[2])]
        try
            b = CI_joint_m_S(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, eff, MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], rel_div_on, f_on, b[3,1]*(1-f_on)/(f_on*Nf_UT), b[2,1]*(1-f_on)/f_on, (1-f_on)*(1+b[2,1])]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], rel_div_on, f_on, b[3,2]*(1-f_on)/(f_on*Nf_UT), b[2,2]*(1-f_on)/f_on, (1-f_on)*(1+b[2,2])]   
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., rel_div_on, f_on, 0., 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, rel_div_on, f_on, Inf, Inf, Inf]  
        end  
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        eff_conv = convert_eff(eff)
        LL_UT = log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, p[1], 1/fit_m[1], eff_conv[1])
        LL_S = log_likelihood_m_S(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], q0_S_off, q_S_off, q0_S_on, q_S_on)
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, fit_m[1], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, p[2], f_on, rel_div_on, fit_m[2], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
	else
		est_res.status .= "failed"
    end   
    return est_res, msel_res
end
# Fraction of on-cells fixed, rel. division rate of on-cells inferred
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Float64, rel_div_on::Bool, fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    start_time = time()
    est_res = est_res_joint_het(cond_S)
    msel_res = msel_res_joint("Heterogeneous", cond_S)
    N_ratio = Nf_S/Nf_UT
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2]) 
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    eff_conv = convert_eff(eff)
    if eff[1] == eff[2]
        eff_conv = (eff_conv, eff_conv)
    end
    # 3 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, rel. division rate on-cells                              
    num_para = 3
    LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], f_on, para[3], q0_UT, q_UT, q0_S_off, q_S_off, 1/fit_m[2], eff_conv[2])
    res = Optim.optimize(LL, [m, S, 1.], iterations=10^4)                                    
    if Optim.converged(res) == true
        est_res.status[5] = "inferred"
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[3], f_on, p[2]*p[1]*(1-f_on)/(f_on*Nf_UT), p[2]*(1-f_on)/f_on, (1-f_on)*(1+p[2])]
        try
            b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], f_on, p[3], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, 1/fit_m[2], eff_conv[2], MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[3,1], f_on, b[4,1]*(1-f_on)/(f_on*Nf_UT), b[2,1]*(1-f_on)/f_on, (1-f_on)*(1+b[2,1])]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[3,2], f_on, b[4,2]*(1-f_on)/(f_on*Nf_UT), b[2,2]*(1-f_on)/f_on, (1-f_on)*(1+b[2,2])] 
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., f_on, 0., 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, f_on, Inf, Inf, Inf]
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m_fitm(mc_counts_UT, mc_max_UT, p[1], 1/fit_m[1], eff_conv[1])
        LL_S = log_likelihood_m_S_div_f(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], f_on, p[3], q0_S_off, q_S_off, 1/fit_m[2], eff_conv[2])
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, fit_m[1], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, p[2], f_on, p[3], fit_m[2], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
	else
		est_res.status .= "failed"
    end
    return est_res, msel_res
end
# Fraction of on-cells not given; inferred if the rel. division rate of on-cells is non-zero
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Bool, rel_div_on::Float64=0., fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    start_time = time()
    est_res = est_res_joint_het(cond_S)[1:5,:]
    msel_res = msel_res_joint("Heterogeneous (zero division rate on-cells)", cond_S)
    N_ratio = Nf_S/Nf_UT
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 1000)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    # For zero rel. division rate on-cells -> Fraction of on-cells cannot be inferred
    if rel_div_on == 0.
        # 2 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio  
        num_para = 2
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
            msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
            LL_UT = log_likelihood_m(mc_counts_UT, mc_max_UT, p[1], q0_UT, q_UT)
            LL_S = log_likelihood_m_S(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], q0_S_off, q_S_off, q0_S_on, q_S_on)
            msel_res.LL = [-MLL, LL_UT, LL_S]
            LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, fit_m[1], eff[1])
            LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, p[2], f_on, rel_div_on, fit_m[2], eff[2])
            msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
            msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
            msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
            end_time = time()
            msel_res.calc_time[1] = end_time - start_time
        else
            est_res.status .= "failed"
        end
    # For non-zero rel. division rate on-cells -> Fraction of on-cells inferred
    else
        # Add the fraction of on-cells and parameters calculated from it
        push!(est_res, ["Fraction on-cells", cond_S, "inferred", 0., 0., Inf])
        push!(est_res, ["Mutation rate on-cells", cond_S, "calc. from 1,4&6", 0., 0., Inf])
        push!(est_res, ["Rel. mutation rate on-cells", cond_S, "calc. from 4&6", 0., 0., Inf])
        push!(est_res, ["Fold change mean mutation rate", cond_S, "calc. from 4&6", 0., 0., Inf])
        msel_res.model .= "Heterogeneous"
        # Calculate the initial value for optimisation
        f_on = initial_f(mc_S, N_ratio, Nf_S, m, S, rel_div_on)
        eff_conv = convert_eff(eff)
        if eff[1] == eff[2]
            eff_conv = (eff_conv, eff_conv)
        end
        # 3 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, fraction of on-cells                              
        num_para = 3
        LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], para[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, 1/fit_m[2], eff_conv[2])
        res = Optim.optimize(LL, [m, S, f_on], iterations=10^4) 
        if Optim.converged(res) == true
            p = Optim.minimizer(res)
            MLL = Optim.minimum(res)
            est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], rel_div_on, p[3], p[2]*p[1]*(1-p[3])/(p[3]*Nf_UT), p[2]*(1-p[3])/p[3], (1-p[3])*(1+p[2])]
            try
                b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], p[3], rel_div_on, q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, false, 1/fit_m[2], eff_conv[2], MLL)
                est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], rel_div_on, b[3,1], b[4,1]/Nf_UT, b[5,1], b[6,1]]
                est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], rel_div_on, b[3,2], b[4,2]/Nf_UT, b[5,2], b[6,2]] 
            catch
                est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., rel_div_on, 0., 0., 0., 0.]
                est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, rel_div_on, 1., Inf, Inf, Inf]

            end
            msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
            LL_UT = log_likelihood_m(mc_counts_UT, mc_max_UT, p[1], q0_UT, q_UT)
            LL_S = log_likelihood_m_S_div_f(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], p[3], rel_div_on, q0_S_off, q_S_off, 1/fit_m[2], eff_conv[2])
            msel_res.LL = [-MLL, LL_UT, LL_S]
            LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, fit_m[1], eff[1])
            LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, p[2], p[3], rel_div_on, fit_m[2], eff[2])
            msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
            msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
            msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
            end_time = time()
            msel_res.calc_time[1] = end_time - start_time
        else
            est_res.status .= "failed"
            end   
        end
    return est_res, msel_res 
end
# Relative division rate on-cells and fraction of on-cells not given -> both inferred
function estimu_het(mc_UT::Vector{Int}, Nf_UT, mc_S::Vector{Int}, Nf_S, eff::Vector{<:Number}, f_on::Bool, rel_div_on::Bool, fit_m::Vector{Float64}=[1., 1.]; cond_S="S")
    start_time = time()
    est_res = est_res_joint_het(cond_S)
    msel_res = msel_res_joint("Heterogeneous", cond_S)
    N_ratio = Nf_S/Nf_UT
    mc_max_UT, mc_counts_UT, num_c_UT = extract_mc(mc_UT)
    mc_max_S, mc_counts_S, num_c_S = extract_mc(mc_S)
    m = max(1.,median(mc_UT))
    S = initial_S(mc_S, m*N_ratio, 10^3)
    f_on = initial_f(mc_S, N_ratio, Nf_S, m, S, 0.)
    q0_UT, q_UT = coeffs(mc_max_UT, 1/fit_m[1], eff[1])
    q0_S_off, q_S_off = coeffs(mc_max_S, 1/fit_m[2], eff[2])
    q0_S_on = -eff[2]
    q_S_on = [eff[2]; zeros(Float64, mc_max_S-1)]
    eff_conv = convert_eff(eff)
    if eff[1] == eff[2]
        eff_conv = (eff_conv, eff_conv)
    end
    # 4 inference parameters: Number of mutations in off-cells under permissive cond., mutation-supply ratio, relative division rate on-cells, fraction of on-cells
    num_para = 4
    LL(para) = -log_likelihood_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, para[1], para[2], para[4], para[3], q0_UT, q_UT, q0_S_off, q_S_off, fit_m[2], eff_conv[2])
    res = Optim.optimize(LL, [m, S, 1., f_on], iterations=10^4)                                                         
    if Optim.converged(res) == true
        est_res.status[5:6] = ["inferred", "inferred"]
        p = Optim.minimizer(res)
        MLL = Optim.minimum(res)
        est_res.MLE = [p[1]/Nf_UT, fit_m[1], fit_m[2], p[2], p[3], p[4], p[2]*p[1]*(1-p[4])/(p[4]*Nf_UT), p[2]*(1-p[4])/p[4], (1-p[4])*(1+p[2])]                                                          
        try
            b = CI_joint_m_S_div_f(mc_counts_UT, mc_max_UT, mc_counts_S, mc_max_S, N_ratio, p[1], p[2], p[4], p[3], q0_UT, q_UT, q0_S_off, q_S_off, q0_S_on, q_S_on, true, 1/fit_m[2], eff_conv[2], MLL)
            est_res.lower_bound = [b[1,1]/Nf_UT, fit_m[1], fit_m[2], b[2,1], b[3,1], b[4,1], b[5,1]/Nf_UT, b[6,1], b[7,1]]
            est_res.upper_bound = [b[1,2]/Nf_UT, fit_m[1], fit_m[2], b[2,2], b[3,2], b[4,2], b[5,2]/Nf_UT, b[6,2], b[7,2]] 
        catch
            est_res.lower_bound = [0., fit_m[1], fit_m[2], 0., 0., 0., 0., 0., 0.]
            est_res.upper_bound = [Inf, fit_m[1], fit_m[2], Inf, Inf, 1., Inf, Inf, Inf]
        end
        msel_res[1,3:5] = selection_crit(num_para, num_c_UT+num_c_S, MLL)
        LL_UT = log_likelihood_m(mc_counts_UT, mc_max_UT, p[1], q0_UT, q_UT)
        LL_S = log_likelihood_m_S_div_f(mc_counts_S, mc_max_S, p[1]*N_ratio, p[2], p[4], p[3], q0_S_off, q_S_off, 1/fit_m[2], eff_conv[2])
        msel_res.LL = [-MLL, LL_UT, LL_S]
        LLs_UT, mc_cutoff_UT, p_cutoff_UT = LL_dist(num_c_UT, Nf_UT, p[1]/Nf_UT, fit_m[1], eff[1])
        LLs_S, mc_cutoff_S, p_cutoff_S = LL_dist(num_c_S, Nf_S, p[1]/Nf_UT, p[2], p[4], p[3], fit_m[2], eff[2])
        msel_res.p_value = 1 .- [ecdf(LLs_UT.+LLs_S)(MLL), ecdf(LLs_UT)(-LL_UT), ecdf(LLs_S)(-LL_S)]
        msel_res.cutoff = [-1, mc_cutoff_UT, mc_cutoff_S]
        msel_res.tail_prob = [NaN, p_cutoff_UT, p_cutoff_S] 
        end_time = time()
        msel_res.calc_time[1] = end_time - start_time
    else
        est_res.status .= "failed"
    end
    return est_res, msel_res
end

function convert_eff(eff::Vector{<:Number})
    # Convert the efficiency vector to a single value or a tuple
    if eff[1] == eff[2] == 1    # Case 1: both efficiencies are equal to 1
        eff_conv = false
    elseif eff[1] == eff[2]
        if eff[1] < 0.5         # Case 2: both efficiencies are equal but less than 0.5
            eff_conv = (eff[1], true)
        else                    # Case 3: both efficiencies are equal and greater than or equal to 0.5 
            eff_conv = eff[1]
        end
    else                        # Case 4: efficiencies are different
        eff_UT, eff_S = eff
        if eff[1] < 0.5
            eff_UT = (eff[1], true)
        end
        if eff[2] < 0.5
            eff_S = (eff[2], true)
        end
        eff_conv = (eff_UT, eff_S)
    end
    return eff_conv
end

extract_mc(mc::Vector{Int}) = maximum(mc), counts(mc, 0:maximum(mc)), length(mc)
selection_crit(num_para, num_c, MLL) = [
    2*(num_para + MLL),                             # AIC
    2*(num_para*(num_c)/(num_c-num_para-1) + MLL),  # AIC corrected (for small sample size) 
    log(num_c)*num_para + 2*MLL                     # BIC
]