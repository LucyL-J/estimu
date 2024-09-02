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
- Install julia separately (instructions and be found here: https://julialang.org/downloads/) and use `julia_setup(JULIA_HOME = "the folder that contains julia binary")`
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
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff=1, fit_m=1., 
       f_on=FALSE, rel_div_on=0., mod="selection", criterion="AIC")
```
This function takes the following inputs.
### Mandatory input
- `mc_UT`: A numeric vector of mutant counts observed in the untreated control condition.
- `Nf_UT`: The average final population size observed in the untreated control condition. Can also be passed as a numeric vector of observed final population sizes, but note that only the average will be used in the inference.
- `mc_S`: A numeric vector of mutant counts observed in the stressed condition (for example, treatment with antimicrobials).
- `Nf_S`: The average final population size observed in the stressed condition. Can also be passed as a numeric vector of observed final population sizes, but see above.
- `plateff`: The plating efficiency, i.e. the fraction of the cultures that is plated onto the selective plates. By default set to $1$ and must be between $>0$ and $\leq 1$. Can be passed as a single value for both conditions, or as a numeric vector with two entries (untreated, stressed) for the untreated and stressed conditions, respectively.

### Optional input
- `fit_m`: The differential fitness of mutant cells defined as the ratio of growth rates of mutant over non-mutant cells ($<1$ if mutants grow slower than non-mutants and $>1$ if they grow faster). By default set to $1$ and must be positive. \
If an estimate of the differential mutant fitness is known from separate experiments, it can be set as a fixed value in the inference. Can be passed as a single value for both conditions, or as a numeric vector with two entries (untreated, stressed) for the untreated and stressed conditions, respectively. \
The differential mutant fitness can  also be inferred. Either, as a joint parameter, i.e. the same value under untreated and stressful conditions, via `fit_m = FALSE`. Or, as two separate inference parameters via `fit_m = c(FALSE, FALSE)`.
- `f_on`: The fraction of the subpopulation with high expression of the stress response ('on-cells'). Applies to the stressed condition only and is assumed to be negligibly small in the untreated control condition. No value assigned by default. \
If an estimate of the fraction of on-cells is known from separate experiments, it can be set as a fixed value in the inference. 
Alternatively, the fraction of on-cells can be inferred via `f_on = FALSE`.
- `rel_div_on`: The relative division rate of on-cells defined as the ratio of division rate of cells with low expression of the stress response ('off-cell') over on-cells. Applies to the stressed condition only and should be $<1$ (by default set to $=0$ based on the SOS response in *E. coli*) because our model is build on the assumption that the expression of the stress response decreases the division rate of a cell. \
If an estimate of the relative division rate of on-cells is known from separate experiments, it can be set as a fixed value in the inference. 
Alternatively, the relative division rate of on-cells can be inferred via `rel_div_on = FALSE`.
- `mod`: The estimation model, see more below.
- `criterion`: If model selection is performed, the criterion that is used to select between homogeneous- and heterogeneous-response model. By default set to the AIC, but can be changed to `criterion = "BIC"`. 

### Estimation under a specific model
The input `mod` can be used to specify a model version that is used in the inference.
1. `null`: Estimation under the null model, i.e. assuming that there is no change in mutation rate between the untreated and the stressed condition. \
The differential fitness of mutants can either be set to fixed value(s) $\neq 1$, or inferred as a joint parameter.
2. `homogeneous`: Estimation under the homogeneous-response model, i.e. assuming that all cells of a population have a comparable level of stress-response expression. \
The mutation rates, and optionally differential mutant fitness, under untreated and stressed conditions are estimated and compared. Values given for `f_on` and `rel_div_on` are ignored. 
3. `heterogeneous`: Estimation under the heterogeneous-response model, i.e. assuming that only a subpopulation highly expressed the stress response which leads to a change in mutation and division rate.\
The differential fitness of mutants can be set to fixed value(s) $\neq 1$, but cannot be inferred. \
For the fraction and relative division rate of on-cells, there are several possibilities: \
**(i)** The fraction of on-cells is known from separate experiments and set to a fixed value in the inference. In this case, the relative division rate of on-cells can either also be set to a fixed value (if known from separate experiments), or inferred.\
If the fraction of on-cells is not known, then, \
**(ii)** if the relative division rate of on-cells is assumed to be $=0$ (the default value), the fraction of on-cells cannot be inferred and is not taken into account in the inference. \
**(iii)** if the relative division rate of on-cells is known from separate experiments and it is $>0$, the fraction of on-cells should be set as an inference parameter. \
**(iv)** if the relative division rate of on-cells is not known, it can either be set to $=0$ (see above), or it can be inferred together with the fraction of on-cells. Note however, that inference of both the fraction and the relative division rate of on-cells is not precise due to model unidentifiability.
4. `standard`: Estimation of the mutation rate, and optionally the differential mutant fitness, for the mutant counts and final population size given as the first two inputs. The second set of mutant counts and final population size is ignored and no comparison between conditions performed. This corresponds to standard mutation rate estimation methods. This mode is also used if the second set of mutant counts or final population size is missing.

### Model selection
When the input is set to `mod = "selection"`, parameters are inferred under null, homogeneous- and heterogeneous-response models the best model/model version is selected.
- If value(s) for the differential mutant fitness are passed, only model versions with fixed differential mutant fitness are compared.
- If no value for the differential mutant fitness is passed, model versions with differential mutant fitness set $=1$ and versions with it inferred as a joint (and as two separate parameters in the case of the homogeneous-response model) are compared.
- If no value for the fraction of on-cells is passed, the heterogeneous-response model version with relative division rate of on-cells set $=0$, and the version with both the relative division and the fraction of on-cells inferred are considered in the model selection (and compared with null and homogeneous-response model)
- If a value for the fraction of on-cells is passed, it is fixed in the inference under the heterogeneous-response. The model version with the relative division rate of on-cells set $=0$, and the version with the relative division rate of on-cells inferred are considered in the model selection.
- Model selection between nested models/model versions is done using likelihood ratio tests (LRTs); and using the AIC/BIC for models that are not nested.

### Estimation output
For each model/model version used in the inference, two dataframes are returned. 

The first dataframe contains information about the model parameters and their estimated values, and has the following columns:
1. `parameter`: Description of the model parameter.
2. `condition`: Condition, to which the parameter applies. Can be a combination of untreated (UT) and stressed (S). For example, the fold change in mutation rate has the condition `S/UT` for 'stressed compared to untreated'.
3. `status`: Whether the parameter is set to a fixed input value, (jointly) inferred, or calculated from other paramters.
4. `MLE`: Maximum-likelihood estimate of the parameter (except for fixed parameters, in which case the input value is given here).
5. `lower_bound`: Lower bound on the parameter estimate, calculated using a profile likelihood approach (except for fixed parameters, in which case the input value is repeated here).
6. `upper_bound`: Uppre bound on the parameter estimate, calculated using a profile likelihood approach (except for fixed parameters, in which case the input value is repeated here).

The second dataframe has the following columns:
1. `model`: The specific model version used in the inference.
2. `selection_result`: If model selection is performed, its result is given here; the options are
- `ns` for 'no support', i.e. there is no significant support for this model version compared to a less complex model version (as determined by an LRT).
- `best hom./het.` for best homogeneous/heterogeneous-response model, i.e. this model version is selection within all homogeneous/heterogeneous-reponse model versions (as determined by an LRT).
- `selected` for the model version selected overall (as determined by the AIC/BIC)
3. `LL` for loglikelihood of the model fit.
4. `AIC` for the AIC value of the model fit.
5. `BIC` for the BIC value of the model fit

## Examples
Input the observed mutant counts, average final population sizes and the plating efficiency
```
mc_UT <- c(80, 0, 9, 3, 11, 0, 0, 1, 0, 3, 5, 1, 1, 2, 0, 0, 1)
Nf_UT <- 1E09
mc_S <- c(10, 2, 4, 0, 0, 0, 0, 1, 1, 2, 1, 1, 0, 2, 1, 0, 1)
Nf_S <- 1.6E08
plateff <- 1
```
Note that the plating efficiency can also be set to different values for the untreated and stressed conditions, respectively. For example,
```
plateff <- c(1, 0.5)
```
Note also that more than one value can be input for the final population sizes. However, only the average values are used in the inference. As an example,
```
Nf_UT <- c(1.6E09, 9E08, 9.5E08)
Nf_S <- c(1.5E08, 1E08, 2E08)
```
### Homogeneous-response model
When estimating under the homogeneous-response model, the differential mutant fitness is set $=1$ by default, and executing the estimation function 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="homogeneous")
```
returns the following
```
[1] "Model used for inference: Homogeneous"
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 9.964485e-10 5.301000e-10 1.659978e-09
2            Mutant fitness        UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
3             Mutation rate         S       inferred 5.049334e-09 2.663162e-09 8.495935e-09
4            Mutant fitness         S   set to input 1.000000e+00 1.000000e+00 1.000000e+00
5      Ratio mutant fitness      S/UT   set to input 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT calc. from 1&3 5.067331e+00 2.234407e+00 1.147827e+01

[[2]]
        model selection_result        LL      AIC      BIC
1 Homogeneous                - -71.63627 147.2725 143.2725
```
The mutation rate under the untreated condition is estimated as $1.0\cdot 10^{-9}$ between $5.3\cdot 10^{-10}$ and $1.7\cdot 10^{-9}$, and the mutation rate under the stressed condition as $5.0\cdot 10^{-9}$ between $2.7\cdot 10^{-9}$ and $8.5\cdot 10^{-9}$. From these mutation rates the fold change in mutation rate is calculated as $5.1$ between $2.2$ and $14.8$. This means that the mutation rate is estimated to be around $5.1$-fold higher in the stressed as in the untreated control condition. \
Note that the lower and upper bounds of the fold change are calculated using a profile-likelihood approach and not by simply dividing the lower/upper bounds of the mutation rate estimates S/UT.

The differential mutant fitness can also be set to a different fixed value (for both the untreated and the stressed conditions) in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="homogeneous", fit_m = 0.8)
```
or as two fixed values, different for the untreated and the stressed conditions, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="homogeneous", fit_m = c(0.8, 0.6))
```

It is also possible to infer the differential mutant fitness, either as a joint inference parameter via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="homogeneous", fit_m = FALSE)
```
which returns the following
```
[1] "Model used for inference: Homogeneous (constr. mutant fitness)"
[[1]]
                  parameter condition           status          MLE  lower_bound  upper_bound
1             Mutation rate        UT         inferred 1.050952e-09 5.573958e-10 1.753498e-09
2            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
3             Mutation rate         S         inferred 5.211676e-09 2.752077e-09 8.753446e-09
4            Mutant fitness      UT+S jointly inferred 7.702790e-01 4.627993e-01 1.310349e+00
5      Ratio mutant fitness                    constr. 1.000000e+00 1.000000e+00 1.000000e+00
6 Fold change mutation rate      S/UT   calc. from 1&3 4.959006e+00 2.205198e+00 1.114531e+01

[[2]]
                                 model selection_result        LL      AIC      BIC
1 Homogeneous (constr. mutant fitness)                - -71.15921 148.3184 152.8975
```
In this case, the differential mutant fitness is constrained to be equal under untreated and stressed conditions, and is estimated as $0.77$ between $0.46$ and $1.31$.

To infer the differential mutant fitness separately under untreated and stressed conditions, execute
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="homogeneous", fit_m = c(FALSE, FALSE))
```
and see the output
```
[1] "Model used for inference: Homogeneous (unconstr. mutant fitness)"
[[1]]
                  parameter condition         status          MLE  lower_bound  upper_bound
1             Mutation rate        UT       inferred 9.867440e-10 5.161329e-10 1.669117e-09
2            Mutant fitness        UT       inferred 1.050062e+00 5.553820e-01 2.125351e+00
3             Mutation rate         S       inferred 5.676540e-09 3.008222e-09 9.498186e-09
4            Mutant fitness         S       inferred 4.010447e-01 1.366183e-01 1.022935e+00
5      Ratio mutant fitness      S/UT calc. from 2&4 3.819249e-01 1.093237e-01 1.146350e+00
6 Fold change mutation rate      S/UT calc. from 1&3 5.752799e+00 2.523235e+00 1.310361e+01

[[2]]
                                   model selection_result        LL      AIC      BIC
1 Homogeneous (unconstr. mutant fitness)                - -69.79779 147.5956 146.6483
```
Now, two different values for the differential mutant fitness are inferred (one for the untreated and one for the stressed condition). Moreover, in addition to the fold change in mutaiton rate, the ratio of mutant fitnesses is calculated; in this case as $0.38$ between $0.11$ and $1.15$.

### Heterogeneous-response model
When estimating under the heterogeneous-response model, by default it is assumed that the fraction of on-cells is not known and the relative division rate of on-cells is set $=0$. Executing the estimation function 
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous")
```
prints the following output
```
[1] "Warning: Fraction of the response-on cells cannot be inferred for zero division rate of response-on cells."
[1] "Model used for inference: Heterogeneous (zero division rate on-cells)"
[[1]]
                    parameter condition           status          MLE  lower_bound  upper_bound
1     Mutation rate off-cells      UT+S jointly inferred 1.032803e-09 5.678150e-10 1.681530e-09
2              Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3              Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4       Mutation-supply ratio         S         inferred 4.880225e+00 1.690687e+00 1.179641e+01
5 Rel. division rate on-cells         S     set to input 0.000000e+00 0.000000e+00 0.000000e+00

[[2]]
                                        model selection_result        LL      AIC      BIC
1 Heterogeneous (zero division rate on-cells)                - -69.44947 142.8989 145.9517
```
In this case, the mutation rate of off-cells and the mutation-supply ratio are inferred. The latter describes how many mutations arise in the subpopulation of on-cells compared to the rest of the population, and is estimated as $4.9$ between $1.6$ and $11.8$ here. This means that around $4.9$ times more mutations arise in on-cells than in the rest of the population. \
Note that for unknown fraction of on-cells and relative division rate of on-cells set $=0$, it is not possible to infer the explicit mutation rate of on-cells or the specific increase in mutation rate due to the induction of the stress response.

In principle, it is possible to infer the fraction of on-cells if either the relative division rate of on-cells is set to fixed value $\neq 0$, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", f_on = FALSE, rel_div_on = 0.1)
```
Or if both fraction and relative division rate of on-cells are inferred via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", f_on = FALSE, rel_div_on = FALSE)
```
which returns the following
```
[1] "Note: Inference of the fraction and relative division rate of response-on cells together not precise."
[1] "Model used for inference: Heterogeneous"
[[1]]
                       parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate off-cells      UT+S jointly inferred 1.032809e-09 5.678150e-10 1.681530e-09
2                 Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3                 Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4          Mutation-supply ratio         S         inferred 4.879880e+00 1.690687e+00 1.179641e+01
5         Mutation rate on-cells         S calc. from 1,4&6 1.119381e-07 2.160171e-09 8.899343e-09
6              Fraction on-cells         S         inferred 4.308487e-02 0.000000e+00 1.000000e+00
7    Rel. division rate on-cells         S         inferred 1.128962e-09 0.000000e+00 9.181609e-01
8    Rel. mutation rate on-cells         S   calc. from 4&6 1.083822e+02 2.050608e-08 8.423223e+58
9 Fold change mean mutation rate      S/UT   calc. from 4&6 5.626546e+00 1.062967e+00 8.423223e+58

[[2]]
          model selection_result        LL      AIC      BIC
1 Heterogeneous                - -69.44947 146.8989 153.0044
```
However in these cases, the estimation of the fraction and relative division rate of on-cells, and all other parameters calculated from these two is not precise due to model unidentifiability. 

The fraction of on-cells can also be set to a fixed value in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", f_on = 0.05)
```
which results in the following
```
[1] "Model used for inference: Heterogeneous (zero division rate on-cells)"
[[1]]
                       parameter condition           status          MLE  lower_bound  upper_bound
1        Mutation rate off-cells      UT+S jointly inferred 1.032803e-09 5.678150e-10 1.681530e-09
2                 Mutant fitness        UT     set to input 1.000000e+00 1.000000e+00 1.000000e+00
3                 Mutant fitness         S     set to input 1.000000e+00 1.000000e+00 1.000000e+00
4          Mutation-supply ratio         S         inferred 4.880225e+00 1.690687e+00 1.179641e+01
5         Mutation rate on-cells         S calc. from 1,4&6 9.576595e-08 4.129363e-08 1.688635e-07
6              Fraction on-cells         S     set to input 5.000000e-02 5.000000e-02 5.000000e-02
7    Rel. division rate on-cells         S     set to input 0.000000e+00 0.000000e+00 0.000000e+00
8    Rel. mutation rate on-cells         S   calc. from 4&6 9.272428e+01 3.212305e+01 2.241318e+02
9 Fold change mean mutation rate      S/UT   calc. from 4&6 5.586214e+00 2.556153e+00 1.215659e+01

[[2]]
                                        model selection_result        LL      AIC      BIC
1 Heterogeneous (zero division rate on-cells)                - -69.44947 142.8989 145.9517
```
Now, the lower and upper bounds on the estimated mutation rate of on-cells ($9.6\cdot 10^{-8}$ between $4.1\cdot 10^{-8}$ and $1.7\cdot 10^{-7}$), and on the estimated relative mutation rate of on-cells, i.e. the specific increase in mutation rate associated with the induction of the stress response, ($92.3$ between $32.1$ and $224.1$) are much tighter than for unknown fraction of on-cells.

The relative division rate of on-cells can also be set to another value $\neq 0$, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", f_on = 0.05, rel_div_on = 0.1)
```
or it can be inferred via
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", f_on = 0.05, rel_div_on = FALSE)
```

The differential mutant fitness can also be set to a different fixed value (for both the untreated and the stressed conditions) in the inference, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", fit_m = 0.8)
```
or as two fixed values, different for the untreated and the stressed conditions, for example
```
estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff, mod="heterogeneous", fit_m = c(0.8, 0.6))
```
Note however, that the differential mutant fitness cannot be inferred for the heterogeneous-response model.