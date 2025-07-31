using StatsBase

# Explaination of the input: mc_max::Int = -1 --> It is simply a dummy value to indicate that no input is given for mc_max

# Explanation of the first line: global max_mc_cutoff = (mc_max == -1) ? max_mc_cutoff_default : max_mc_cutoff_default*log10(mc_max+1) --> Concise if-else assignment
    # If no input is given for the maximal obseved mutant count mc_max it is set to -1 and the default value defined in line 14 in estimationfunctions.jl is used for the cutoff
    # Else, the value of mc_max is used to scale the default value 
# This way, a larger cutoff is considered for data with larger counts and vice versa
# The choice of log10 is arbitrary, but it ensures the cutoff is not too large for large mc_max values
# Previously, I set the cutoff to a fixed value independent of mc_max, which led to large tail probabilities and made GoF testing less reliable when counts were large
# I also tried using a fixed tail probability (of 0.05) instead, but the results were less good than the current approach

function LL_dist(num_c::Int, N, mu, fit_m, eff, mc_max::Int = -1)
    global max_mc_cutoff = (mc_max == -1) ? max_mc_cutoff_default : max_mc_cutoff_default*log10(mc_max+1)
    mc_rdraws, p_rdraws = r_mudi(R_gof*num_c, N, mu, fit_m, eff)
    LL = Vector{Float64}(undef, R_gof)
    s = randperm(R_gof * num_c) 
    c_indices = reshape(s, num_c, R_gof)
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
    end
    return LL, length(p_rdraws), p_rdraws[end]
end

function LL_dist(num_c::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff, mc_max::Int = -1)
    global max_mc_cutoff = (mc_max == -1) ? max_mc_cutoff_default : max_mc_cutoff_default*log10(mc_max+1)
    mc_rdraws, p_rdraws = r_mudi(R_gof*num_c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    LL = Vector{Float64}(undef, R_gof)
    s = randperm(R_gof * num_c) 
    c_indices = reshape(s, num_c, R_gof)
    for i in eachindex(LL)
        mc_obs = mc_rdraws[c_indices[:, i]]
        mc_counts_obs = counts(mc_obs, 0:maximum(mc_obs))
        p = p_rdraws[1:maximum(mc_obs)+1]
        LL[i] = -sum(mc_counts_obs .* log.(p)) 
    end
    return LL, length(p_rdraws), p_rdraws[end]
end