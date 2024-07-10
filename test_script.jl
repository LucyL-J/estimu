include("estimationfunctions.jl")

mc_UT = [80, 0, 9, 3, 11, 0, 0, 1, 0, 3, 5, 1, 1, 2, 0, 0, 1];
Nf_UT = 10^9;
mc_S = [10, 2, 4, 0, 0, 0, 0, 0, 1, 2, 1, 0, 0, 2, 1, 0, 1];
Nf_S = 1.6*10^8;

for eff_UT in [1, 0.9, 0.1]
    for fit_m in (1., 0.8, false)
        println(eff_UT)
        println(estimu(mc_UT, Nf_UT/eff_UT, eff_UT, fit_m)[1])
    end
end

for eff_UT in [1, 0.9]
    for eff_S in [1, 0.1]
        for fit_m in [[1.,1.], [0.8,0.4]]
            println((eff_UT, eff_S))
            println(estimu_0(mc_UT, Nf_UT/eff_UT, mc_S, Nf_S/eff_S, [eff_UT,eff_S], fit_m)[1])
            for f_on in (0.05, false)
                for rel_div_on in (0., 0.1, false)
                    println((eff_UT, eff_S))
                    println(estimu_het(mc_UT, Nf_UT/eff_UT, mc_S, Nf_S/eff_S, [eff_UT,eff_S], f_on, rel_div_on, fit_m)[1])
                end
            end
        end
        for fit_m in ([1.,1.], [0.8,0.4], false, (false,false))
            println((eff_UT, eff_S))
            println(estimu_hom(mc_UT, Nf_UT/eff_UT, mc_S, Nf_S/eff_S, [eff_UT,eff_S], fit_m)[1])
        end
    end
end
