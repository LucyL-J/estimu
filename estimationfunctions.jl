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
        q0 = q0_coeff(1/fit_m, eff)
        q = q_coeffs(mc_max, 1/fit_m, eff)
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