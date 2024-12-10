using StatsBase
   
function r_mudi(K::Int, N, mu, fit_m, eff)
    R_mudi = zeros(Int, K)
    m_max_guess = Int(round(N*mu*fit_m*eff*10))
    P = p_mudi(m_max_guess, N, mu, fit_m, eff)
    k = 1
    p = P[k]
    R = rand(K)
    sort!(R)
    for j in eachindex(R)
        r = R[j]
        while r > p
        k += 1
        if k > m_max_guess
            m_max_guess *= 10
            P = p_mudi(m_max_guess, N, mu, fit_m, eff)
        end
        p += P[k]
        end
        R_mudi[j] = k-1
    end
    return R_mudi
end

function r_mudi(K::Int, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    R_mudi = zeros(Int, K)
    m_max_guess = Int(round(N*mu_off*((1-f_on)*fit_m + S*f_on)*eff*10))
    P = p_mudi(m_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    k = 1
    p = P[k]
    R = rand(K)
    sort!(R)
    for j in eachindex(R)
        r = R[j]
        while r > p
        k += 1
        if k > m_max_guess
            m_max_guess *= 10
            P = p_mudi(m_max_guess, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
        end
        p += P[k]
        end
        R_mudi[j] = k-1
    end
    return R_mudi
end

function KL_divs(R, c, N, mu, fit_m, eff)
    KL_div = zeros(Float64, R)
    mc_total = r_mudi(R*c, N, mu, fit_m, eff)
    P = p_mudi(maximum(mc_total), N, mu, fit_m, eff)
    s = sample(1:R*c, R*c, replace=false)
    for i in eachindex(KL_div)
        mc_sample = mc_total[1+(i-1)*c:i*c]
        p_sample = counts(mc_sample, 0:maximum(mc_sample))
        p = P[collect(1:maximum(mc_sample)+1)][p_sample .!= 0]
        p_sample = p_sample[p_sample .!= 0] ./ length(mc_sample)
        KL_div[i] = sum(p_sample .* log.(p_sample ./ p))
    end
    return KL_div
end

function KL_divs(R, c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    KL_div = zeros(Float64, R)
    mc_total = r_mudi(R*c, N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    P = p_mudi(maximum(mc_total), N, mu_off, S, f_on, rel_div_on, fit_m, eff)
    s = sample(1:R*c, R*c, replace=false)
    for i in eachindex(KL_div)
        mc_sample = mc_total[1+(i-1)*c:i*c]
        p_sample = counts(mc_sample, 0:maximum(mc_sample))
        p = P[collect(1:maximum(mc_sample)+1)][p_sample .!= 0]
        p_sample = p_sample[p_sample .!= 0] ./ length(mc_sample)
        KL_div[i] = sum(p_sample .* log.(p_sample ./ p))
    end
    return KL_div
end



    