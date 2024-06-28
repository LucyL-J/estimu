using SpecialFunctions

# Pre-calculations independent of the inference parameter m
q(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
function q(K::Int, inv_fit_m)
    q = zeros(Float64, K)
    q[1] = inv_fit_m/(inv_fit_m+1)
    for k = 2:K
        @inbounds q[k] = (k-1)/(k+inv_fit_m) * q[k-1]
    end
    return q
end   
scale_f(f_on, rel_div_on) = (1-f_on)/(1-f_on*(1-rel_div_on))
inverse_fit_on(f_on, rel_div_on) = (1-f_on*(1-rel_div_on))/rel_div_on

# Mutant count distributions (probabilities to obverse 0:K mutant colonies)

# When mutants have zero differential mutant fitness -> Poisson distribution
mudi_0(m, K::Int, f) = [m^k for k = 0:K] ./ f .* exp(-m)

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# When mutants have non-zero differential mutant fitness
# q = rf for differential fitness = 1 and q = f/g else
function mudi(m, K::Int, q0, q)
    p = zeros(Float64, K+1)
    p[1] = exp(-m*q0)
    # Recursive calculation of probabilities
    for k = 1:K
        @views S = sum((1:k) .* q[1:k] .* reverse(p[1:k]))
        if S > 0.
            @inbounds p[k+1] = m*S/k
        end
    end
    return p
end

# Probability density and cumulative distribution functions

# Poisson (fit_m=0), Luria-Dellbrueck (fit_m=1) or Mandelbrot-Koch else
function p_mudi(K::Int, N, mu, fit_m=1.)
    if fit_m == 0.
        f = factorials(K)
        p = mudi_0(mu*N, K, f)
    else
        if fit_m == 1.
            q = reduced_fs(K)
        else 
            q = factorials(K-1) ./ gammas(K, 1/fit_m)
        end
        p = mudi(mu*N, K, q)
    end
    return p
end
pmf_mudi(k::Int, N, mu, fit_m=1.) = p_mudi(k, N, mu, fit_m)[k+1]
cdf_mudi(k::Int, N, mu, fit_m=1.) = sum(p_mudi(k, N, mu, fit_m))

# With a subpopulation of on-cells
function ps_mudi(K::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m=1.)
    if rel_div_on == 0.
        # on-cells have zero division rate
        f = factorials(K)
        p_on = mudi_0(mu_on*N*f_on/(1-f_on), K, f)
        if fit_m == 1.
            rf = reduced_fs(K)
            p_off = mudi(mu_off*N, K, rf)
        else
            g = gammas(K, 1/fit_m)
            p_off = mudi(mu_off*N, K, f[1:end-1]./g)
        end
    else
        f = factorials(K-1)
        # Re-scale the effective population size (because on-cells now also contribute to growth)
        N *= scale_f(f_on, rel_div_on)
        # Calculate differential fitness of on-cells (compared to total population growth)
        ifit = inverse_fit_on(f_on, rel_div_on)/fit_m
        g_on = gammas(K, ifit)
        p_on = mudi(mu_on*N*f_on/(1-f_on), K, f./g_on)
        if fit_m == 1.
            rf = reduced_fs(K)
            p_off = mudi(mu_off*N, K, rf)
        else
            g_off = gammas(K, 1/fit_m)
            p_off = mudi(mu_off*N, K, f./g_off)
        end
    end
    return p_off, p_on
end
function pmf_mudi(k::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m=1.)
    p_off, p_on = ps_mudi(k, N, mu_off, mu_on, f_on, rel_div_on, fit_m)
    return sum(p_off .* reverse(p_on))    # Total mutant count = sum of contributions from on-cells and the rest of the population
end
function cdf_mudi(k::Int, N, mu_off, mu_on, f_on, rel_div_on, fit_m=1.)
    p_off, p_on = ps_mudi(k, N, mu_off, mu_on, f_on, rel_div_on, fit_m)
    return sum([sum(p_off[1:i] .* reverse(p_on[1:i])) for i = 1:k+1])
end
