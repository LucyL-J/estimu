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
