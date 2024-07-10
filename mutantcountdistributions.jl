using Distributions, HypergeometricFunctions

# Pre-calculations independent of the inference parameter m

# Plating efficiency = 1 without differential mutant fitness (= 1)
q_coeffs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
# Plating efficiency = 1 with differential mutant fitness
function q_coeffs(K::Int, inv_fit_m)
    q = zeros(Float64, K)
    q[1] = inv_fit_m/(inv_fit_m+1)
    for k = 2:K
        @inbounds q[k] = (k-1)/(k+inv_fit_m) * q[k-1]
    end
    return q
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
function mudi(K::Int, m, q0, q)
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
mudi(K::Int, m) = pdf(Poisson(m), 0:K)
mudi(K::Int, m, eff) = pdf(Poisson(m*eff), 0:K)
# With subpopulation of on-cells
function mudi(K::Int, m_off, q0_off, q_off, m_on, q0_on, q_on)
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
        p = mudi(K, mu*N, eff)
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
        p = mudi(K, mu*N, q0, q)
    end
    return p
end
pmf_mudi(k::Int, N, mu, fit_m, eff=1.) = p_mudi(max(2,k), N, mu, fit_m, eff)[k+1]
cdf_mudi(k::Int, N, mu, fit_m, eff=1.) = sum(p_mudi(max(2,k), N, mu, fit_m, eff)[1:k+1])

# With a subpopulation of on-cells
function p_mudi(K::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m, eff=1.)
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
    return mudi(K, N*mu_off, q0_off, q_off, N*mu_on*f_on/(1-f_on), q0_on, q_on)
end
pmf_mudi(k::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m, eff=1.) = p_mudi(max(2,k), N, mu_off, mu_on, f_on, rel_div_on, fit_m, eff)[k+1]
cdf_mudi(k::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m, eff=1.) = sum(p_mudi(max(2,k), N, mu_off, mu_on, f_on, rel_div_on, fit_m, eff)[1:k+1])
