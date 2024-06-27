using SpecialFunctions

# Pre-calculations independent of the inference parameter m
factorials(K::Int) = [factorial(big(k)) for k in 0:K-1]
reduced_fs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
gammas(K::Int, inv_fit_m) = [gamma(1+inv_fit_m + k) for k in 1:K] ./ gamma(1+inv_fit_m) ./ inv_fit_m
