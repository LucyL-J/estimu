using SpecialFunctions

# Pre-calculations independent of the inference parameter m
factorials(K::Int) = Float64.([factorial(big(k)) for k in 0:K])
reduced_fs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
gammas(K::Int, inv_fit_m) = [gamma(1+inv_fit_m + k) for k in 1:K] ./ gamma(1+inv_fit_m) ./ inv_fit_m

# Mutant count distributions (probabilities to obverse 0:K mutant colonies)

# When mutants have zero differential mutant fitness -> Poisson distribution
mudi_0(m, K::Int, f) = [m^k for k = 0:K] ./ f .* exp(-m)

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011

# When mutants have non-zero differential mutant fitness
# q = rf for differential fitness = 1 and q = f/g else
function mudi(m, K::Int, q)
    p = zeros(Float64, K+1)
    p[1] = exp(-m)
    # Recursive calculation of probabilities
    q .*= m
    for k = 1:K
        @views S = sum(reverse(1:k) .* reverse(q[1:k]) .* p[1:k])
        if S > 0.
            @inbounds p[k+1] = S/k
        end
    end
    return p
end

# Probability density and cumulative distribution SpecialFunctions

function pdf_mudi(k::Int, mu, N, fit_m=1.)
    if fit_m == 0.
        f = factorials(k)
        p = mudi_0(mu*N, k, f)
    else
        if fit_m == 1.
            q = reduced_fs(k)
        else 
            q = factorials(k-1) ./ gammas(k, 1/fit_m)
        end
        p = mudi(mu*N, k, q)
    end
    return p[k+1]
end