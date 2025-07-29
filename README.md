# estimu
ESTimating Increases in Mutation rates

## Download and installation
1. Install R, and RStudio (recommended), by following the instructions here: https://rstudio-education.github.io/hopr/starting.html

2. Launch R/RStudio
3. Install the package JuliaCall via executing
```
install.packages("JuliaCall")
``` 
4. Install julia
- This can be done through the package JuliaCall via executing
```
library(JuliaCall)
julia_setup(installJulia = TRUE)
```
- Alternatively, install julia separately, by following the instructions here: https://julialang.org/downloads/, and then executing
```
library(JuliaCall)
julia_setup()
```
5. Download this repository. Click on the dropdown menu `<> Code` on the top right above the files and 
- either, `Download ZIP` (if you do not have any experience with git)
- or, copy the url to clipboard and clone the repository with `git clone <url>` from the terminal

## Setup of estimu
1. Open the file **setup.R** from the downloaded/cloned repository in R/RStudio
2. Set the working directory to the repository folder
- In RStudio, rightclick on the file **setup.R** and select `Set Working Directory`
- In the console, execute `setwd("the folder that contains setup.R file")`
3. Run the setup.R file via executing
```
source("setup.R")
```

### Troubleshooting
1. Error in julia_setup() : julia is not found
- Install julia separately (instructions and be found here: https://julialang.org/downloads/) and use `julia_setup(JULIA_HOME = "the folder that contains julia binary")`. You can run Sys.BINDIR in julia to get the full path to the directory containing the julia binary.
- More information on the package JuliaCall can be found here: https://cran.r-project.org/web/packages/JuliaCall/readme/README.html
2. If there is still an error, try again after installing R version 4.4.1 and julia version 1.10.4 (this tool was written on these specific versions)
3. It might be necessary to activate the environment and install the packages within julia to avoid errors because of dependencies. For this, install and launch julia in a terminal. Navigate to the folder of the repository via `cd("the folder that contains the file Project.toml")`, then execute,
```
using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
```

## Running estimu
Estimation of increases in mutation rates is done via the function
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff=1, 
       fit_m=1., f_on=FALSE, rel_div_on=0., mod="all")
```
This function takes the following inputs.
### Mandatory input
- `mc_UT`: A numeric vector of mutant counts observed in the untreated control condition.
- `Nf_UT`: The average final population size observed in the untreated control condition. Can also be passed as a numeric vector of observed final population sizes, but note that only the average will be used in the inference.
- `mc_S`: A numeric vector of mutant counts observed in the stressed condition (for example, treatment with antimicrobials).
- `Nf_S`: The average final population size observed in the stressed condition. Can also be passed as a numeric vector of observed final population sizes, but see above.
- `plateff`: The plating efficiency, i.e. the fraction of the cultures that is plated onto the selective plates. By default set to $1$ and must be between $0$ and $1$. Can be passed as a single value for both conditions, or as a numeric vector with two entries (untreated, stressed) for the untreated and stressed conditions, respectively.

### Optional input
- `mod`: The estimation model, see more in the next section.
- `fit_m`: The differential fitness of mutant cells defined as the ratio of growth rates of mutant over non-mutant cells ($<1$ if mutants grow slower than non-mutants and $>1$ if they grow faster). By default set to $1$ and must be positive. \
If an estimate of the differential mutant fitness is known from separate experiments, it can be set as a fixed value in the inference. Can be passed as a single value for both conditions, or as a numeric vector with two entries (untreated, stressed) for the untreated and stressed conditions, respectively. \
For homogeneously expressed stress-responses, the differential mutant fitness can  also be inferred. Either, as a joint parameter, i.e. the same value under untreated and stressful conditions, via `fit_m = FALSE`. Or as two separate inference parameters via `fit_m = c(FALSE, FALSE)`.
#### Input specific to heterogeneously expressed stress responses:
- `f_on`: The fraction of the subpopulation with high expression of the stress response ('on-cells'). Applies to the stressed condition only and is assumed to be negligibly small in the untreated control condition. No value assigned by default. \
If an estimate of the fraction of on-cells is known from separate experiments, it can be set as a fixed value in the inference. 
Alternatively, the fraction of on-cells can be inferred via `f_on = FALSE`.
- `rel_div_on`: The relative division rate of on-cells defined as the ratio of division rate of cells with low expression of the stress response ('off-cell') over on-cells. Applies to the stressed condition only and should be $<1$ (by default set to $=0$ based on the SOS response in *E. coli*) because our model is build on the assumption that the expression of the stress response decreases the division rate of a cell. \
If an estimate of the relative division rate of on-cells is known from separate experiments, it can be set as a fixed value in the inference. 
Alternatively, the relative division rate of on-cells can be inferred via `rel_div_on = FALSE`.

### Estimation under a specific model
The input `mod` can be used to specify a model version that is used in the inference.
1. `null`: Estimation under the null model, i.e. assuming that there is no change in mutation rate between the untreated and the stressed condition. \
The differential fitness of mutants can either be set to fixed value(s), or inferred as one joint or two separate inference parameter.
2. `hom`: Estimation under the homogeneous-response model, i.e. assuming that all cells of a population have a comparable level of stress-response expression. \
The mutation rates, and optionally differential mutant fitness, under untreated and stressed conditions are estimated separately and compared. Values given for `f_on` and `rel_div_on` are ignored. 
3. `het`: Estimation under the heterogeneous-response model, i.e. assuming that only a subpopulation highly expressed the stress response which leads to a change in mutation and division rate.\
The differential fitness of mutants can be set to fixed value(s) $\neq 1$, but cannot be inferred. \
For the fraction and relative division rate of on-cells, there are several possibilities: \
**(i)** The fraction of on-cells is known from separate experiments and set to a fixed value in the inference. In this case, the relative division rate of on-cells can either also be set to a fixed value (if known from separate experiments), or inferred.\
If the fraction of on-cells is not known, then, \
**(ii)** if the relative division rate of on-cells is assumed to be $=0$ (the default value), the fraction of on-cells cannot be inferred and is not taken into account in the inference. \
**(iii)** if the relative division rate of on-cells is known from separate experiments and it is $>0$, the fraction of on-cells should be set as an inference parameter. \
**(iv)** if the relative division rate of on-cells is not known, it can either be set to $=0$ (see above), or it can be inferred together with the fraction of on-cells. Note however, that inference of both the fraction and the relative division rate of on-cells together is not precise due to model unidentifiability.
4. `standard`: Estimation of the mutation rate, and optionally the differential mutant fitness, for a single set of mutant counts and final population size (the first two inputs are used). The second set of mutant counts and final population size is ignored and no comparison between conditions performed. This corresponds to standard mutation rate estimation methods. This mode is also used if the second set of mutant counts or final population size is missing.

### Comparing all models
When the input is set to `mod = "all"`, parameters are inferred under all null, homogeneous- and heterogeneous-response models using fixed values specified in the input. The AIC, AIC<sub>c</sub>, BIC, log-likelihood and goodness-of-fit p-values can be used for model comparison and selection.
- If value(s) for the differential mutant fitness are passed, only model versions with fixed differential mutant fitness are compared.
- If no value for the differential mutant fitness is passed, model versions with differential mutant fitness set $=1$ and versions with it inferred as a joint and as two separate inference parameters are compared.
- If no value for the relative division rate of on-cells is passed, the heterogeneous-response model version with relative division rate of on-cells set $=0$, and the version with both the relative division and the fraction of on-cells inferred are considered.
- If a value for the relative division of on-cells is passed, the heterogeneous-response model version with fixed relative division rate and inferred fraction of on-cells, and the version with the relative division rate of on-cells set $=0$, are considered.

### Estimation output
For each model/model version used in the inference, two dataframes are returned. The dataframes can be accessed via saving the output and using the double square bracket syntax.
```
res <- estimu(...)
res[[1]]
res[[2]]
```

The first dataframe contains information about the model parameters and their estimated values, and has the following columns:
1. `parameter`: Description of the model parameter.
2. `condition`: Condition, to which the parameter applies. Can be a combination of untreated (`UT`) and stressed (`S`). For example, the fold change in mutation rate has the condition `S/UT` for 'stressed compared to untreated'.
3. `status`: Whether the parameter is set to a fixed input value, (jointly) inferred, or calculated from other paramters.
4. `MLE`: Maximum-likelihood estimate of the parameter (except for fixed parameters, in which case the input value is given here).
5. `lower_bound`: Lower bound on the parameter estimate, calculated using a profile likelihood approach (except for fixed parameters, in which case the input value is repeated here).
6. `upper_bound`: Uppre bound on the parameter estimate, calculated using a profile likelihood approach (except for fixed parameters, in which case the input value is repeated here).

The second dataframe has the following columns:
1. `model`: The specific model version used in the inference.
2. `condition`: Condition, to which the subsequent values belong: Untreated (`UT`), stressed (`S`) or joint (`UT+S`). In cases with independent estimation for the untreated and the stressed condition, a fourth 'condition' (`UT->S`) is reported. It gives the results of a goodness-of-fit test of whether the data input for the stressed condition could have been generated under the parameters estimated for the untreated condition.
3. `AIC`: The AIC value of the model fit(s), defined as AIC $:= -2 \ln L + 2p$ with $\ln L$ being the loglikelihood and $p$ the number of inference parameters.
4. `AIC_corr`: The AIC<sub>c</sub> (corrected for small sample sizes) value of the model fit(s). The AIC<sub>c</sub> is defined as $:= -2 \ln L + 2p \frac{n}{n-p-1}$ with $n$ the sample size, i.e. the total number of parallel cultures in the untreated and stressed conditions together.
5. `BIC`: The BIC value of the model fit(s), defined as BIC $:= -2 \ln L + p \ln n$. \
Note that the AIC, AIC<sub>c</sub> and BIC are not defined for the individual conditions (`UT` and `S`) if there is at least one joint inference parameter.
6. `LL`: The loglikelihoods of the model fits.
7. `p_value`: The p-value of the goodness-of(-model)-fit tests. A value $<0.05$ indicates a poor fit of the model to the specific condition. For the special case of (`UT->S`), a value $<0.05$ indicates that the data input for the stressed condition is significantly different from the data input for the untreated condition. Conversely, a value $\geq 0.05$ hints towards no significant difference between conditions.
8. `cutoff`: The maximal mutant count at which the simulation-based goodness-of-fit test is cutoff. Only applicable to individual and not joint conditions (otherwise, a value of `-1` is returned).
9. `tail_prob`: The probability of observing a mutant count larger then the above-described cutoff. A large value here indicates computational problems in the goodness-of-fit test and the p-values might not be reliable. 
10. `calc_time`: The time required for parameter estimation and goodness-of-fit tests in seconds.

## Examples
Input the observed mutant counts, average final population sizes and the plating efficiency
```
mc_UT <- c(80, 0, 9, 3, 11, 0, 0, 1, 0, 3, 5, 1, 1, 2, 0, 0, 1)
Nf_UT <- 1.15E09
mc_S <- c(10, 2, 4, 0, 0, 0, 0, 1, 1, 2, 1, 1, 0, 2, 1, 0, 1)
Nf_S <- 1.5E08
plateff <- 1
```
Note that the plating efficiency can also be set to different values for the untreated and stressed conditions, respectively. For example,
```
plateff <- c(1, 0.5)
```
Note also that more than one value can be input for the final population sizes. However, only the average values are used in the inference. As an example, the following input will result in the exact same estimates as the input from above because the average final population sizes are the same:
```
Nf_UT <- c(1.6E09, 9E08, 9.5E08)
Nf_S <- c(1.5E08, 1E08, 2E08)
```
### Homogeneous-response model
When estimating under the homogeneous-response model, the differential mutant fitness is set $=1$ by default, and executing the estimation function 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="hom")
```
returns the following
```
[1] "Model used for inference: Homogeneous"
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 8.664769e-10 4.609565e-10 1.443459e-09
2            Mutant fitness        UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
3             Mutation rate         S       inferred 5.385957e-09 2.840707e-09 9.062331e-09
4            Mutant fitness         S   set to input 1.000000e+00 1.000000e+00 1.000000e+00
5      Ratio mutant fitness      S/UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT calc. from 1&3 6.215926e+00 3.278456e+00 1.168431e+01

[[2]]
        model condition       AIC  AIC_corr       BIC        LL p_value cutoff    tail_prob calc_time
1 Homogeneous      UT+S 147.27254 147.65964 150.32526 -71.63627  0.6782     -1          NaN 0.1752172
2 Homogeneous        UT  84.96144  85.22811  85.79465 -41.48072  0.4254   1909 2.754722e-07 0.0479970
3 Homogeneous         S  62.31110  62.57776  63.14431 -30.15555  0.7963   1042 7.517840e-07 0.1126490
4 Homogeneous     UT->S       NaN       NaN       NaN -30.43262  0.0060   1042 1.200002e-07 0.1752172
```
The mutation rate under the untreated condition is estimated as $8.7\cdot 10^{-10}$ between $4.6\cdot 10^{-10}$ and $1.4\cdot 10^{-9}$, and the mutation rate under the stressed condition as $5.4\cdot 10^{-9}$ between $2.8\cdot 10^{-9}$ and $9.1\cdot 10^{-9}$. From these mutation rates the fold change in mutation rate is calculated as $6.2$ between $3.3$ and $11.7$. This means that the mutation rate is estimated to be around $6.2$-fold higher in the stressed as in the untreated control condition. \
The model fits the data adequately for untreated, stressed and joint conditions (respective p-values $>0.05$), and the stressed condition is significantly different from the untreated condition (respective `UT->S` p-value $<0.05$).

The differential mutant fitness can also be set to a different fixed value (for both the untreated and the stressed conditions) in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="hom", fit_m = 0.8)
```
or as two fixed values, different for the untreated and the stressed conditions, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="hom", fit_m = c(0.8, 0.6))
```

It is also possible to infer the differential mutant fitness, either as a joint inference parameter via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="hom", fit_m = FALSE)
```
which returns the following
```
[1] "Model used for inference: Homogeneous (constr. mutant fitness)"
[[1]]
                  parameter condition           status          MLE  lower_bound  upper_bound
1             Mutation rate        UT         inferred 9.138711e-10 4.846920e-10 1.524781e-09
2            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
3             Mutation rate         S         inferred 5.559122e-09 2.935549e-09 9.337009e-09
4            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
5      Ratio mutant fitness      S/UT          constr. 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT   calc. from 1&3 6.083048e+00 3.212213e+00 1.146939e+01

[[2]]
                                 model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Homogeneous (constr. mutant fitness)      UT+S 148.3184 149.1184 152.8975 -71.15921  0.4512     -1          NaN  1.103412
2 Homogeneous (constr. mutant fitness)        UT      NaN      NaN      NaN -41.90725  0.2275   1909 4.610239e-08       NaN
3 Homogeneous (constr. mutant fitness)         S      NaN      NaN      NaN -29.25195  0.7192   1042 1.474308e-07       NaN
```
In this case, the differential mutant fitness is constrained to be equal under untreated and stressed conditions, and is estimated as $0.77$ between $0.46$ and $1.31$.

To infer the differential mutant fitness separately under untreated and stressed conditions, execute
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="hom", fit_m = c(FALSE, FALSE))
```
and see the output
```
[1] "Model used for inference: Homogeneous (unconstr. mutant fitness)"
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 8.580382e-10 4.488112e-10 1.451406e-09
2            Mutant fitness        UT       inferred 1.050062e+00 5.553820e-01 2.125351e+00
3             Mutation rate         S       inferred 6.054975e-09 3.208770e-09 1.013140e-08
4            Mutant fitness         S       inferred 4.010447e-01 1.366183e-01 1.022935e+00
5      Ratio mutant fitness      S/UT calc. from 2&4 3.819249e-01 1.301051e-01 9.741662e-01
6 Fold change mutation rate      S/UT calc. from 1&3 7.056766e+00 3.739659e+00 1.349114e+01

[[2]]
                                   model condition       AIC  AIC_corr       BIC        LL p_value cutoff    tail_prob  calc_time
1 Homogeneous (unconstr. mutant fitness)      UT+S 147.59557 148.97488 153.70102 -69.79779  0.4715     -1          NaN 0.33614516
2 Homogeneous (unconstr. mutant fitness)        UT  86.94050  87.79764  88.60692 -41.47025  0.4582   1909 3.654037e-07 0.22175908
3 Homogeneous (unconstr. mutant fitness)         S  60.65508  61.51222  62.32150 -28.32754  0.4638    135 2.798150e-07 0.01877403
4 Homogeneous (unconstr. mutant fitness)     UT->S       NaN       NaN       NaN -30.62331  0.0077   1042 1.545957e-07 0.33614516
```
Now, two different values for the differential mutant fitness are inferred (one for the untreated and one for the stressed condition). Moreover, in addition to the fold change in mutaiton rate, the ratio of mutant fitnesses is calculated; in this case as $0.38$ between $0.13$ and $0.97$, implying that mutant are less fit than non-mutants and that this effect is more pronounced in the stressed than the untreated condition.

### Heterogeneous-response model
When estimating under the heterogeneous-response model, by default it is assumed that the fraction of on-cells is not known and the relative division rate of on-cells is set $=0$. Executing the estimation function 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het")
```
gives the following output
```
[1] "Warning: Fraction of the response-on cells cannot be inferred for zero division rate of response-on cells."
[1] "Model used for inference: Heterogeneous (zero division rate on-cells)"
[[1]]
                    parameter condition           status          MLE  lower_bound  upper_bound
1     Mutation rate off-cells      UT+S jointly inferred 9.072963e-10 4.984535e-10 1.478130e-09
2              Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3              Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4       Mutation-supply ratio         S         inferred 6.157943e+00 2.252714e+00 1.457909e+01
5 Rel. division rate on-cells         S     set to input 0.000000e+00 0.000000e+00 0.000000e+00

[[2]]
                                        model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Heterogeneous (zero division rate on-cells)      UT+S 143.1065 143.4936 146.1592 -69.55324  0.3827     -1          NaN  0.324368
2 Heterogeneous (zero division rate on-cells)        UT      NaN      NaN      NaN -41.49372  0.4695   1909 2.885438e-07       NaN
3 Heterogeneous (zero division rate on-cells)         S      NaN      NaN      NaN -28.05953  0.2812   1042 1.258652e-07       NaN
```
In this case, the mutation rate of off-cells and the mutation-supply ratio are inferred. The latter describes how many mutations arise in the subpopulation of on-cells compared to the rest of the population, and is estimated as $6.2$ between $2.6$ and $14.8$ here. This means that around $6.2$ times more mutations arise in the subpopulation of on-cells than in the rest of the population. \
Note that for unknown fraction of on-cells and relative division rate of on-cells set $=0$, it is not possible to infer the explicit mutation rate of on-cells or the specific increase in mutation rate due to the induction of the stress response.

In principle, it is possible to infer the fraction of on-cells if either the relative division rate of on-cells is set to fixed value $\neq 0$, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", f_on = FALSE, rel_div_on = 0.1)
```
Or if both fraction and relative division rate of on-cells are inferred via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", f_on = FALSE, rel_div_on = FALSE)
```
which returns the following
```
[1] "Note: Inference of the fraction and relative division rate of response-on cells together not precise."
[1] "Model used for inference: Heterogeneous"
[[1]]
                       parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate off-cells      UT+S jointly inferred 9.073221e-10 4.984534e-10 1.478130e-09
2                 Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3                 Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4          Mutation-supply ratio         S         inferred 6.158557e+00 2.252714e+00 1.457909e+01
5    Rel. division rate on-cells         S         inferred 5.577405e-09 0.000000e+00 9.521684e-01
6              Fraction on-cells         S         inferred 6.199004e-02 0.000000e+00 1.000000e+00
7         Mutation rate on-cells         S calc. from 1,4&6 8.455243e-08 1.316071e-18          Inf
8    Rel. mutation rate on-cells         S   calc. from 4&6 9.318898e+01 1.450500e-09          Inf
9 Fold change mean mutation rate      S/UT   calc. from 4&6 6.714798e+00 1.686026e-09 1.367533e+01

[[2]]
          model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Heterogeneous      UT+S 147.1065 148.4858 153.2119 -69.55324  0.3766     -1          NaN 0.9534001
2 Heterogeneous        UT      NaN      NaN      NaN -41.49373  0.4719   1909 2.885521e-07       NaN
3 Heterogeneous         S      NaN      NaN      NaN -28.05951  0.2783   1042 1.258688e-07       NaN
```
However in these cases, the estimation of the fraction and relative division rate of on-cells, and all other parameters calculated from these two is not precise due to model unidentifiability. (The 95% confidence intervals span several orders of magnitude, or the entire interval of possible values from $0$ to $1$ in the case of the fraction of on-cells.)

The fraction of on-cells can also be set to a fixed value in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", f_on = 0.05)
```
which results in the following
```
[1] "Model used for inference: Heterogeneous (zero division rate on-cells)"
[[1]]
                       parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate off-cells      UT+S jointly inferred 9.072963e-10 4.984535e-10 1.478130e-09
2                 Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3                 Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4          Mutation-supply ratio         S         inferred 6.157943e+00 2.252714e+00 1.457909e+01
5    Rel. division rate on-cells         S     set to input 0.000000e+00 0.000000e+00 0.000000e+00
6              Fraction on-cells         S     set to input 5.000000e-02 5.000000e-02 5.000000e-02
7         Mutation rate on-cells         S calc. from 1,4&6 1.061545e-07 3.883369e-08 2.513235e-07
8    Rel. mutation rate on-cells         S   calc. from 4&6 1.170009e+02 4.280156e+01 2.770026e+02
9 Fold change mean mutation rate      S/UT   calc. from 4&6 6.800046e+00 3.090078e+00 1.480013e+01

[[2]]
                                        model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Heterogeneous (zero division rate on-cells)      UT+S 143.1065 143.4936 146.1592 -69.55324  0.3698     -1          NaN 0.1353889
2 Heterogeneous (zero division rate on-cells)        UT      NaN      NaN      NaN -41.49372  0.4658   1909 2.885438e-07       NaN
3 Heterogeneous (zero division rate on-cells)         S      NaN      NaN      NaN -28.05953  0.2716   1042 1.258652e-07       NaN
```
Now, the lower and upper bounds on the estimated mutation rate of on-cells ($1.1\cdot 10^{-7}$ between $3.9\cdot 10^{-8}$ and $2.5\cdot 10^{-7}$), and on the estimated relative mutation rate of on-cells, i.e. the specific increase in mutation rate associated with the induction of the stress response, ($117$ between $43$ and $277$) are tighter than for unknown fraction of on-cells.

The relative division rate of on-cells can also be set to another value $> 0$, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", f_on = 0.05, rel_div_on = 0.1)
```
or it can be inferred via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", f_on = 0.05, rel_div_on = FALSE)
```

The differential mutant fitness can also be set to a different fixed value (for both the untreated and the stressed conditions) in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", fit_m = 0.8)
```
or as two fixed values, different for the untreated and the stressed conditions, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="het", fit_m = c(0.8, 0.6))
```
Note however, that the differential mutant fitness cannot be set as an inference parameter for the heterogeneous-response model.

### Model comparison and selection
To compare different null, homogeneous- and heterogeneous-response models set the input `mod="all"`, which estimates parameters under all possible model versions in accordance with fixed input parameters. The AIC, AIC<sub>c</sub> or BIC values can be used to compare the models and select the best (with the lowest respective value). Executing 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all")
```
gives the following output
```
[1] "Estimated parameters under all models are:"
[[1]]
[[1]]
             parameter condition           status          MLE lower_bound  upper_bound
1        Mutation rate      UT+S jointly inferred 1.437369e-09 9.34704e-10 2.082672e-09
2       Mutant fitness        UT     set to input 1.000000e+00 1.00000e+00 1.000000e+00
3       Mutant fitness         S     set to input 1.000000e+00 1.00000e+00 1.000000e+00
4 Ratio mutant fitness      S/UT     set to input 1.000000e+00 1.00000e+00 1.000000e+00

[[2]]
  model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1  Null      UT+S 162.3972 162.5222 163.9236 -80.19862  0.0960     -1          NaN 0.1344459
2  Null        UT      NaN      NaN      NaN -43.36645  0.8755   1909 4.590640e-07       NaN
3  Null         S      NaN      NaN      NaN -36.83217  0.0068   1042 1.992613e-07       NaN

[[2]]
[[1]]
             parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate      UT+S jointly inferred 1.500129e-09 9.715099e-10 2.180045e-09
2       Mutant fitness        UT jointly inferred 7.643513e-01 4.612581e-01 1.290671e+00
3       Mutant fitness         S jointly inferred 7.643513e-01 4.612581e-01 1.290671e+00
4 Ratio mutant fitness      S/UT          constr. 1.000000e+00 1.000000e+00 1.000000e+00

[[2]]
                          model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Null (constr. mutant fitness)      UT+S 163.3607 163.7478 166.4134 -79.68035  0.0248     -1          NaN 0.2645209
2 Null (constr. mutant fitness)        UT      NaN      NaN      NaN -43.77926  0.6152   1909 7.132875e-08       NaN
3 Null (constr. mutant fitness)         S      NaN      NaN      NaN -35.90109  0.0028   1042 3.741706e-08       NaN

[[3]]
[[1]]
             parameter condition              status          MLE  lower_bound  upper_bound
1        Mutation rate      UT+S    jointly inferred 1.481083e-09 9.552621e-10 2.159605e-09
2       Mutant fitness        UT            inferred 8.910953e-01 4.774560e-01 1.750107e+00
3       Mutant fitness         S            inferred 5.727896e-01 2.364955e-01 1.372375e+00
4 Ratio mutant fitness      S/UT calculated from 2&3 6.427928e-01 2.653986e-01 1.540099e+00

[[2]]
                            model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Null (unconstr. mutant fitness)      UT+S 164.7037 165.5037 169.2828 -79.35187  0.0368     -1          NaN 0.8383899
2 Null (unconstr. mutant fitness)        UT      NaN      NaN      NaN -43.55883  0.7916   1909 2.227449e-07       NaN
3 Null (unconstr. mutant fitness)         S      NaN      NaN      NaN -35.79304  0.0008    500 2.422268e-08       NaN

[[4]]
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 8.664769e-10 4.609565e-10 1.443459e-09
2            Mutant fitness        UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
3             Mutation rate         S       inferred 5.385957e-09 2.840707e-09 9.062331e-09
4            Mutant fitness         S   set to input 1.000000e+00 1.000000e+00 1.000000e+00
5      Ratio mutant fitness      S/UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT calc. from 1&3 6.215926e+00 3.278456e+00 1.168431e+01

[[2]]
        model condition       AIC  AIC_corr       BIC        LL p_value cutoff    tail_prob  calc_time
1 Homogeneous      UT+S 147.27254 147.65964 150.32526 -71.63627  0.6690     -1          NaN 0.06332731
2 Homogeneous        UT  84.96144  85.22811  85.79465 -41.48072  0.4280   1909 2.754722e-07 0.02606320
3 Homogeneous         S  62.31110  62.57776  63.14431 -30.15555  0.7899   1042 7.517840e-07 0.02338910
4 Homogeneous     UT->S       NaN       NaN       NaN -30.43262  0.0065   1042 1.200002e-07 0.06332731

[[5]]
[[1]]
                  parameter condition           status          MLE  lower_bound  upper_bound
1             Mutation rate        UT         inferred 9.138711e-10 4.846920e-10 1.524781e-09
2            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
3             Mutation rate         S         inferred 5.559122e-09 2.935549e-09 9.337009e-09
4            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
5      Ratio mutant fitness      S/UT          constr. 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT   calc. from 1&3 6.083048e+00 3.212213e+00 1.146939e+01

[[2]]
                                 model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Homogeneous (constr. mutant fitness)      UT+S 148.3184 149.1184 152.8975 -71.15921  0.4582     -1          NaN 0.2989171
2 Homogeneous (constr. mutant fitness)        UT      NaN      NaN      NaN -41.90725  0.2211   1909 4.610239e-08       NaN
3 Homogeneous (constr. mutant fitness)         S      NaN      NaN      NaN -29.25195  0.7290   1042 1.474308e-07       NaN

[[6]]
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 8.580382e-10 4.488112e-10 1.451406e-09
2            Mutant fitness        UT       inferred 1.050062e+00 5.553820e-01 2.125351e+00
3             Mutation rate         S       inferred 6.054975e-09 3.208770e-09 1.013140e-08
4            Mutant fitness         S       inferred 4.010447e-01 1.366183e-01 1.022935e+00
5      Ratio mutant fitness      S/UT calc. from 2&4 3.819249e-01 1.301051e-01 9.741662e-01
6 Fold change mutation rate      S/UT calc. from 1&3 7.056766e+00 3.739659e+00 1.349114e+01

[[2]]
                                   model condition       AIC  AIC_corr       BIC        LL p_value cutoff    tail_prob  calc_time
1 Homogeneous (unconstr. mutant fitness)      UT+S 147.59557 148.97488 153.70102 -69.79779  0.4651     -1          NaN 0.20814919
2 Homogeneous (unconstr. mutant fitness)        UT  86.94050  87.79764  88.60692 -41.47025  0.4533   1909 3.654037e-07 0.07634401
3 Homogeneous (unconstr. mutant fitness)         S  60.65508  61.51222  62.32150 -28.32754  0.4579    312 1.467656e-08 0.01777506
4 Homogeneous (unconstr. mutant fitness)     UT->S       NaN       NaN       NaN -30.62331  0.0073   1042 1.545957e-07 0.20814919

[[7]]
[[1]]
                    parameter condition           status          MLE  lower_bound  upper_bound
1     Mutation rate off-cells      UT+S jointly inferred 9.072963e-10 4.984535e-10 1.478130e-09
2              Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3              Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4       Mutation-supply ratio         S         inferred 6.157943e+00 2.252714e+00 1.457909e+01
5 Rel. division rate on-cells         S     set to input 0.000000e+00 0.000000e+00 0.000000e+00

[[2]]
                                        model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Heterogeneous (zero division rate on-cells)      UT+S 143.1065 143.4936 146.1592 -69.55324  0.3731     -1          NaN 0.1291051
2 Heterogeneous (zero division rate on-cells)        UT      NaN      NaN      NaN -41.49372  0.4617   1909 2.885438e-07       NaN
3 Heterogeneous (zero division rate on-cells)         S      NaN      NaN      NaN -28.05953  0.2810   1042 1.258652e-07       NaN

[[8]]
[[1]]
                       parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate off-cells      UT+S jointly inferred 9.073221e-10 4.984534e-10 1.478130e-09
2                 Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3                 Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4          Mutation-supply ratio         S         inferred 6.158557e+00 2.252714e+00 1.457909e+01
5    Rel. division rate on-cells         S         inferred 5.577405e-09 0.000000e+00 9.521684e-01
6              Fraction on-cells         S         inferred 6.199004e-02 0.000000e+00 1.000000e+00
7         Mutation rate on-cells         S calc. from 1,4&6 8.455243e-08 1.316071e-18          Inf
8    Rel. mutation rate on-cells         S   calc. from 4&6 9.318898e+01 1.450500e-09          Inf
9 Fold change mean mutation rate      S/UT   calc. from 4&6 6.714798e+00 1.686026e-09 1.367533e+01

[[2]]
          model condition      AIC AIC_corr      BIC        LL p_value cutoff    tail_prob calc_time
1 Heterogeneous      UT+S 147.1065 148.4858 153.2119 -69.55324  0.3810     -1          NaN 0.4283109
2 Heterogeneous        UT      NaN      NaN      NaN -41.49373  0.4690   1909 2.885521e-07       NaN
3 Heterogeneous         S      NaN      NaN      NaN -28.05951  0.2834   1042 1.258688e-07       NaN
```

All models expect the null models are adequate fits to the data, indicating that the mutation rate is indeed different between the untreated and stressed condition.\
Since the sample size, i.e. total number of parallel cultures, is $<40$ we recommend using the corrected AIC<sub>c</sub> for model selection. The model with the lowest AIC<sub>c</sub> is the heterogeneous-response model with zero relative division rate of on-cells. The model with the second lowest AIC<sub>c</sub> is the homogeneous-response model without differential mutant fitness. The difference in AIC<sub>c</sub> is roughly $\Delta$ AIC $=4.2$ suggesting considerable support for a heterogeneously expressed stress response.

In the example above, neither differential mutant fitness nor fraction of on-cells were passed as input. Therefore, the model versions with default fixed values and with inferred parameters were compared. 

For the null model, these model versions are
1. Differential mutant fitness set $=1$ for untreated and stressed conditions
2. Differential mutant fitness inferred as a joint parameter (constrained to be equal under untreated and stressed conditions)
3. Differential mutant fitness inferred separately for untreated and stressed conditions (unconstrained)

For the homogeneous-response model, the model versions are
1. Differential mutant fitness set $=1$ for untreated and stressed conditions
2. Differential mutant fitness inferred as a joint parameter (constrained to be equal under untreated and stressed conditions)
3. Differential mutant fitness inferred separately for untreated and stressed conditions (unconstrained)

For the heterogeneous-response model, the model versions are
1. Relative division rate of on-cells set $=0$ and fraction of on-cells not taken into account
2. Relative division rate and fraction of on-cells both inferred

It is also possible to specify the differential mutant fitness and/or the relative division rate of on-cells when estimating under all models, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all", fit_m = 0.8)
```
or
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all", fit_m = c(0.8, 0.6))
```
or 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all", rel_div_on = 0.1)
```
or any combination of the above
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all", fit_m = 0.8, rel_div_on = 0.1)
```
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="all", fit_m = c(0.8, 0.6), rel_div_on = 0.1)
```