# Calculations of the mutant count distributions are based on 
# Łazowski, K. (2023). Efficient, robust, and versatile fluctuation data analysis using MLE MUtation Rate calculator (mlemur). Mutation Research - Fundamental and Molecular Mechanisms of Mutagenesis, 826(April). https://doi.org/10.1016/j.mrfmmm.2023.111816

using Distributions, HypergeometricFunctions, Random

max_mc_cutoff = 1000 # Maximum number of mutants to be considered in the mutant count distribution

# Pre-calculations independent of the inference parameter m

# Plating efficiency = 1 without differential mutant fitness (= 1)
q_coeffs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
function q_coeffs!(q, k)
    push!(q, 1/(k*(k+1)))
end
# Plating efficiency = 1 with differential mutant fitness
function q_coeffs(K::Int, inv_fit_m)
    q = zeros(Float64, K)
    q[1] = inv_fit_m/(inv_fit_m+1)
    for k = 2:K
        @inbounds q[k] = (k-1)/(k+inv_fit_m) * q[k-1]
    end
    return q
end   
function q_coeffs!(q, k, inv_fit_m)
    push!(q, (k-1)/(k+inv_fit_m) * q[k-1])
end
# Partial plating with efficiency < 1
q0_coeff(inv_fit_m, eff) = -1 + inv_fit_m*(1-eff)/(inv_fit_m+1) * pFq((1,1), (inv_fit_m+2,), 1-eff)
function beta_f(K::Int, inv_fit_m)
    B = zeros(Float64, K)
    B[1] = inv_fit_m/(inv_fit_m+1)
    for k = 2:K
        @inbounds B[k] = (k-1)/(k+inv_fit_m) * B[k-1]
    end
    return B
end
# For efficiency > 0.5, calculate hypergeometric functions reversely from K:1
function q_coeffs(K::Int, inv_fit_m, eff::Float64)
    F = zeros(Float64, K)
    F[K] = pFq((inv_fit_m, inv_fit_m+1), (inv_fit_m+1+K, ), 1 - eff)
    F[K-1] = pFq((inv_fit_m, inv_fit_m+1), (inv_fit_m+K, ), 1 - eff)
    for k = K-1:-1:2
        f = inv_fit_m + k
        @inbounds F[k-1] = k*(k+1)*(1-eff)/(eff*f*(f+1)) * F[k+1] + (f - 2*k*(1-eff))/(eff*f) * F[k]
    end
    B = beta_f(K, inv_fit_m)
    return @. eff^inv_fit_m * B * F
end
function q_coeffs!(q, k, b, inv_fit_m, eff)
    inv_fit_m_big = BigFloat(inv_fit_m)
    eff_big = BigFloat(eff)
    push!(q, Float64(eff^inv_fit_m * b * pFq((inv_fit_m_big, inv_fit_m_big+1), (inv_fit_m_big+1+k, ), 1 - eff_big)))
end
# For efficiency < 0.5, calculate hypergeometric functions from 1:K
function q_coeffs(K::Int, inv_fit_m, eff::Tuple{Float64,Bool})
    eff = eff[1]
    F = zeros(Float64, K)
    F[1] = pFq((inv_fit_m, inv_fit_m+1), (inv_fit_m+2, ), 1 - eff)
    F[2] = pFq((inv_fit_m, inv_fit_m+1), (inv_fit_m+3, ), 1 - eff)
    for k = 2:K-1
        f = inv_fit_m + k
        @inbounds F[k+1] = (f+1)/(k*(k+1)*(1-eff)) * (eff*f*F[k-1] - (f-2*k*(1-eff))*F[k]) 
    end
    B = beta_f(K, inv_fit_m)
    return @. eff^inv_fit_m * B * F
end

function coeffs(K::Int, inv_fit_m, eff)
    if eff == 1.
        q0 = -1
        if inv_fit_m == 1.
            q = q_coeffs(K)
        else
            q = q_coeffs(K, inv_fit_m)
        end
    else
        q0 = q0_coeff(inv_fit_m, eff)
        if eff < 0.5
            q = q_coeffs(K, inv_fit_m, (eff, true))
        else
            q = q_coeffs(K, inv_fit_m, eff)
        end
    end
    return q0, q
end

# Factor to rescale effective population size (due to subpopulation of on-cells contributing to population growth)
scale_f(f_on, rel_div_on) = (1-f_on)/(1-f_on*(1-rel_div_on))
# Fitness of on-cells compared to population growth
inverse_fit_on(f_on, rel_div_on) = (1-f_on*(1-rel_div_on))/rel_div_on

# Mutant count distributions (probabilities to obverse 0:K mutant colonies)
function mudi_K(K::Int, m, q0, q)
    p = zeros(Float64, K+1)
    p[1] = exp(m*q0)
    # Recursive calculation of probabilities
    for k = 1:K
        @views S = sum((1:k) .* q[1:k] .* reverse(p[1:k]))
        @inbounds p[k+1] = m*S/k
    end
    return max.(p, 0.)
end
# Zero differential mutant fitness
mudi_K(K::Int, m) = pdf(Poisson(m), 0:K)
mudi_K(K::Int, m, eff) = pdf(Poisson(m*eff), 0:K)
# With subpopulation of on-cells
function mudi_K(K::Int, m_off, q0_off, q_off, m_on, q0_on, q_on)
    p = zeros(Float64, K+1)
    p[1] = exp(m_off*q0_off + m_on*q0_on)
    for k = 1:K
        @views S = sum((1:k) .* (m_off.*q_off[1:k] .+ m_on.*q_on[1:k]) .* reverse(p[1:k]))
        @inbounds p[k+1] = S/k
    end
    return max.(p, 0.)
end

# Probability mass and cumulative distribution functions
# Poisson (fit_m=0), Luria-Dellbrueck (fit_m=1) or Mandelbrot-Koch else
function p_mudi(K::Int, N, mu, fit_m, eff=1.)
    if fit_m == 0.
        p = mudi_K(K, mu*N, eff)
    else
        if eff == 1.
            q0 = -1
            q = q_coeffs(K, 1/fit_m)
        else
            q0 = q0_coeff(1/fit_m, eff)
            if eff < 0.5
                q = q_coeffs(K, 1/fit_m, (eff, true))
            else
                q = q_coeffs(K, 1/fit_m, eff)
            end
        end
        p = mudi_K(K, mu*N, q0, q)
    end
    return p
end
pmf_mudi(k::Int, N, mu, fit_m, eff=1.) = p_mudi(max(2,k), N, mu, fit_m, eff)[k+1]
cdf_mudi(k::Int, N, mu, fit_m, eff=1.) = sum(p_mudi(max(2,k), N, mu, fit_m, eff)[1:k+1])

# With a subpopulation of on-cells
function p_mudi(K::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff=1.)
    if eff == 1.
        q0_off = -1
        q_off = q_coeffs(K, 1/fit_m)
    else
        q0_off = q0_coeff(1/fit_m, eff)
        if eff < 0.5
            q_off = q_coeffs(K, 1/fit_m, (eff, true))
        else
            q_off = q_coeffs(K, 1/fit_m, eff)
        end
    end
    if rel_div_on == 0.
        # on-cells have zero division rate
        q0_on = -eff
        q_on = zeros(Float64, K)
        q_on[1] = eff            
    else
        # Re-scale the effective population size (because on-cells now also contribute to growth)
        N *= scale_f(f_on, rel_div_on)
        # Calculate differential fitness of on-cells (compared to total population growth)
        ifit = inverse_fit_on(f_on, rel_div_on)/fit_m
        if eff == 1.
            q0_on = -1
            q_on = q_coeffs(K, ifit)
        else
            q0_on = q0_coeff(ifit, eff)
            if eff < 0.5
                q_on = q_coeffs(K, ifit, (eff, true))
            else
                q_on = q_coeffs(K, ifit, eff)
            end
        end
    end
    return mudi_K(K, N*mu_off, q0_off, q_off, N*mu_off*S, q0_on, q_on)
end
pmf_mudi(k::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff=1.) = p_mudi(max(2,k), N, mu_off, S, f_on, rel_div_on, fit_m, eff)[k+1]
cdf_mudi(k::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff=1.) = sum(p_mudi(max(2,k), N, mu_off, S, f_on, rel_div_on, fit_m, eff)[1:k+1])

# Mutant count distributions until given threshold of the cdf

# Adding the probability of observing k mutants to the vector of probabilities 0:k-1
# For a homogeneous population
function add_p!(p, q, m, k)
    @views S = sum((1:k) .* q[1:k] .* reverse(p[1:k]))
    push!(p, isnan(m*S/k) ? 0 : max(m*S/k, 0.))
end
# For a heterogeneous population
function add_p!(p, q_off, q_on, m_off, m_on, k)
    @views S = sum((1:k) .* (m_off.*q_off[1:k] .+ m_on.*q_on[1:k]) .* reverse(p[1:k]))
    push!(p, isnan(S/k) ? 0 : max(S/k, 0.))
end

function tail_prob(p, max_mc_cutoff)    # "Tail" probability = 1 - cumulative probability
    if p[end] == 0
        p[end] = 1 - sum(p)
    elseif length(p) == max_mc_cutoff   # "Jackpot" of observing large number of mutants
        push!(p, 1-sum(p)) 
    end
    return p
end

function mudi_threshold(p_threshold, m) # Homogeneous population
    p_cumulative = exp(-m)
    q = Vector{Float64}(undef, 0)
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q, k)
        add_p!(p, q, m, k)
        k += 1
        p_cumulative += p[k]
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold(p_threshold, m, inv_fit_m) # With diff. mutant fitness
    p_cumulative = exp(-m)
    q = [inv_fit_m/(inv_fit_m+1)]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        add_p!(p, q, m, k)
        k += 1
        p_cumulative += p[k]
        q_coeffs!(q, k, inv_fit_m)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold(p_threshold, m, inv_fit_m, eff) # Plating efficiency < 1
    p_cumulative = exp(m*(-1 + inv_fit_m*(1-eff)/(inv_fit_m+1) * pFq((1,1), (inv_fit_m+2,), 1-eff)))
    b = inv_fit_m/(inv_fit_m+1)
    q = Vector{Float64}(undef, 0)
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q, k, b, inv_fit_m, eff)
        add_p!(p, q, m, k)
        k += 1
        b *= (k-1)/(k+inv_fit_m)
        p_cumulative += p[k]
    end
    return tail_prob(p, max_mc_cutoff)
end

function mudi_threshold_het_0(p_threshold, m_off, m_on) # Heterogeneous population, zero division rate of response-on subpopulation
    p_cumulative = exp(-(m_off + m_on))
    q_off = Vector{Float64}(undef, 0)
    q_on = [1.]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q_off, k)
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        p_cumulative += p[k]
        push!(q_on, 0.)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold_het_0(p_threshold, m_off, m_on, inv_fit_m) # With diff. mutant fitness in response-off subpopulation
    p_cumulative = exp(-(m_off + m_on))
    q_off = [inv_fit_m/(inv_fit_m+1)]
    q_on = [1.]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        p_cumulative += p[k]
        push!(q_on, 0.)
        q_coeffs!(q_off, k, inv_fit_m)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold_het_0(p_threshold, m_off, m_on, inv_fit_m, eff) # Plating efficiency < 1
    p_cumulative = exp(m_off*(-1 + inv_fit_m*(1-eff)/(inv_fit_m+1) * pFq((1,1), (inv_fit_m+2,), 1-eff)) - eff)
    b = inv_fit_m/(inv_fit_m+1)
    q_off = Vector{Float64}(undef, 0)
    q_on = [eff]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q_off, k, b, inv_fit_m, eff)
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        b *= (k-1)/(k+inv_fit_m)
        p_cumulative += p[k]
        push!(q_on, 0.)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold_het(p_threshold, m_off, m_on, ifit) # None-zero relative division rate of response-on subpopulation
    p_cumulative = exp(-(m_off + m_on))
    q_off = Vector{Float64}(undef, 0)
    q_on = [ifit/(ifit+1)]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q_off, k)
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        p_cumulative += p[k]
        q_coeffs!(q_on, k, ifit)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold_het(p_threshold, m_off, m_on, inv_fit_m, ifit) # Diff. mutant fitness in response-off subpopulation
    p_cumulative = exp(-(m_off + m_on))
    q_off = [inv_fit_m/(inv_fit_m+1)]
    q_on = [ifit/(ifit+1)]
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        p_cumulative += p[k]
        q_coeffs!(q_off, k, inv_fit_m)
        q_coeffs!(q_on, k, ifit)
    end
    return tail_prob(p, max_mc_cutoff)
end
function mudi_threshold_het(p_threshold, m_off, m_on, inv_fit_m, ifit, eff) # Plating efficiency < 1
    p_cumulative = exp(m_off*(-1 + inv_fit_m*(1-eff)/(inv_fit_m+1) * pFq((1,1), (inv_fit_m+2,), 1-eff)) + m_on*(-1 + ifit*(1-eff)/(ifit+1) * pFq((1,1), (ifit+2,), 1-eff)))
    b_off = inv_fit_m/(inv_fit_m+1)
    b_on = ifit/(ifit+1)
    q_off = Vector{Float64}(undef, 0)
    q_on = Vector{Float64}(undef, 0)
    p = [p_cumulative]
    k = 1
    while p_threshold > p_cumulative && k < max_mc_cutoff && p[k] > 0
        q_coeffs!(q_off, k, b_off, inv_fit_m, eff)
        q_coeffs!(q_on, k, b_on, ifit, eff)
        add_p!(p, q_off, q_on, m_off, m_on, k)
        k += 1
        b_off *= (k-1)/(k+inv_fit_m)
        b_on *= (k-1)/(k+ifit)
        p_cumulative += p[k]
    end
    return tail_prob(p, max_mc_cutoff)
end

function r_mudi(K::Int, N, mu, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    k = 1
    uni_draws = sort(rand(K))
    p_threshold = uni_draws[end]
    if eff == 1.
        if fit_m == 1.
            p_draws = mudi_threshold(p_threshold, N*mu)
        else
            p_draws = mudi_threshold(p_threshold, N*mu, 1/fit_m)
        end
    else
        p_draws = mudi_threshold(p_threshold, N*mu, 1/fit_m, eff)
    end
    cumulative_p = p_draws[1]
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p && k <= max_mc_cutoff 
            k += 1
            cumulative_p += p_draws[k]
        end
        random_draws[j] = k-1
    end
    return random_draws, p_draws
end

function r_mudi(K::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    random_draws = Vector{Int}(undef, K)
    k = 1
    uni_draws = sort(rand(K))
    p_threshold = uni_draws[K]
    if rel_div_on == 0.
        if eff == 1.
            if fit_m == 1.
                p_draws = mudi_threshold_het_0(p_threshold, N*mu_off, N*mu_off*S)
            else
                p_draws = mudi_threshold_het_0(p_threshold, N*mu_off, N*mu_off*S, 1/fit_m)
            end
        else
            p_draws = mudi_threshold_het_0(p_threshold, N*mu_off, N*mu_off*S, 1/fit_m, eff)
        end
    else
        N *= scale_f(f_on, rel_div_on)
        ifit = inverse_fit_on(f_on, rel_div_on)/fit_m
        if eff == 1.
            if fit_m == 1.
                p_draws = mudi_threshold_het(p_threshold, N*mu_off, N*mu_off*S, ifit)
            else
                p_draws = mudi_threshold_het(p_threshold, N*mu_off, N*mu_off*S, 1/fit_m, ifit)
            end
        else
            p_draws = mudi_threshold_het(p_threshold, N*mu_off, N*mu_off*S, 1/fit_m, ifit, eff)
        end
    end
    cumulative_p = p_draws[1]
    for j in eachindex(uni_draws)
        r = uni_draws[j]
        while r > cumulative_p && k <= max_mc_cutoff
            k += 1
            cumulative_p += p_draws[k]
        end
        random_draws[j] = k-1
    end
    return random_draws, p_draws[1:k]
end