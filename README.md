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