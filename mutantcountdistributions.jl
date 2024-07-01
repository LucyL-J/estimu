q_coeffs(K::Int) = 1 ./ [k*(k+1) for k in 1:K]
function q_coeffs(K::Int, inv_fit_m)
q0_coeff(inv_fit_m, eff) = -1 + inv_fit_m*(1-eff)/(inv_fit_m+1) * pFq((1,1), (inv_fit_m+2,), 1-eff)
function q_coeffs(K::Int, inv_fit_m, eff)
# For efficiency < 0.5, calculate hypergeometric functions from 1:K
function q_coeffs(K::Int, inv_fit_m, eff, small_eff::Bool)
function mudi(K::Int, m, q0, q)
# Zero differential mutant fitness
mudi(K::Int, m) = pdf(Poisson(m), 0:K)
mudi(K::Int, m, eff) = pdf(Poisson(m*eff), 0:K)
