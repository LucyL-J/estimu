using SpecialFunctions

# Pre-calculations independent of the inference parameter m
factorials(K::Int) = [factorial(big(k)) for k in 0:K]
reduced_fs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
gammas(K::Int, inv_fit_m) = [gamma(1+inv_fit_m + k) for k in 1:K] ./ gamma(1+inv_fit_m) ./ inv_fit_m

# Mutant count distributions (probabilities to obverse 0:K mutant colonies)

# When mutants have zero differential mutant fitness -> Poisson distribution
mudi_hom_0(m, K::Int, f) = Float64.([m^k for k = 0:K] ./ f .* exp(-m))

# The following recursive formulas for the pdf of the mutant count distribution are based on 
# Keller, P., & Antal, T. (2015). Mutant number distribution in an exponentially growing population. Journal of Statistical Mechanics: Theory and Experiment, 2015(1), P01011. https://doi.org/10.1088/1742-5468/2015/01/P01011