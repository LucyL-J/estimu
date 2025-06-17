library(JuliaCall)
julia_setup()
# Run Sys.BINDIR in julia to get the full path to the directory containing the julia binary
#julia_setup(JULIA_HOME = "/Users/lucylansch-justen/.julia/juliaup/julia-1.11.3+0.aarch64.apple.darwin14/bin")

for (p in c("StatsBase", "Optim", "DataFrames", "Distributions", "HypergeometricFunctions", "Roots")) {
  julia_install_package_if_needed(p)
}
julia_source("estimationfunctions.jl")

check_input <- function(mc, Nf, plateff=1, fit_m=1., rel_div_on=FALSE, f_on=0.1){
  status <- TRUE
  if(is.numeric(mc) && (sum(mc>=0)==length(mc))){
    if(sum(as.integer(mc)!=mc)>0){
      print("Warning: Mutant counts have been rounded to integers for the inference.")
    }
    mc <- as.integer(mc)
  } else {
    status <- FALSE
    print("Error: Mutant counts must be positive numbers.")
  }
  if(is.numeric(Nf) && (sum(Nf>=0)==length(Nf)) && (sum(Nf)>0)){
    if(length(Nf)>1){
      Nf <- mean(Nf)
      print("Note: The mean of the given final population sizes is used in the inference.")
    }
  } else {
    status <- FALSE
    print("Error: Final population sizes must be positive numbers.")
  }
  if(is.numeric(plateff) && (sum(plateff>0)==length(plateff)) && (sum(plateff<=1)==length(plateff))){
    if(length(plateff)==1){
      plateff <- c(plateff, plateff)
    }
    if(length(plateff)>2){
      plateff <- plateff[1:2]
      print("Note: More than two values for plating efficency given; only the first two (untreated & stressful) will be used in the inference.")
    }
  } else {
    status <- FALSE
    print("Error: Plating efficency(s) must be strictly between zero and one.")
  }
  if(is.numeric(fit_m) && (sum(fit_m>=0)==length(fit_m))){
    if(length(fit_m==1)){
      fit_m <- c(fit_m, fit_m)
    } else {
      if(length(fit_m)>2){
        fit_m <- fit_m[1:2]
        print("Note: More than two values for mutant fitness given; only the first two (untreated & stressful) will be used in the inference.")
      } 
    }
  } else {
    if(is.logical(fit_m)){
      if(length(fit_m)>2){
        fit_m <- fit_m[1:2]
        print("Note: More than two values for mutant fitness given; only the first two (untreated & stressful) will be used in the inference.")
      }
    } else {
      status <- FALSE
      print("Error: Mutant fitness must be either given as positive number(s), or set as a joint inference parameter via fit_m=FALSE or as two separate inference parameters via fit_m=c(FALSE,FALSE)")
    } 
  }
  if((is.numeric(rel_div_on) && (rel_div_on >= 0)) || is.logical(rel_div_on)){
  } else {
    status <- FALSE
    print("Error: Relative division rate of response-on cells has to be positive (by default = 0), or set as an inference parameter via rel_div_on=FALSE")
  }
  if(is.logical(f_on)){
    if((is.numeric(rel_div_on))&&(rel_div_on==0)){
      print("Warning: Fraction of the response-on cells cannot be inferred for zero division rate of response-on cells.")
    } else {
      print("Note: Inference of the fraction and relative division rate of response-on cells together not precise.") 
    }
  } else {
    if(is.numeric(f_on) && f_on > 0 && f_on < 100){
      if(f_on >= 1){
        f_on <- f_on[1]/100
        print("Note: Assuming that the fraction of the response-on cells was given in percent.")
      } 
    } else {
      status <- FALSE
      print("Error: Fraction of response-on cells has to be strictly between zero and one.")
    }
  }
  return(list(status, mc, as.numeric(Nf), plateff, fit_m, f_on, rel_div_on))
}

estimu <- function(mc_UT, Nf_UT, mc_S, Nf_S, plateff=1, fit_m=1., f_on=FALSE, rel_div_on=0., mod="selection", criterion="AIC"){
  res <- "Warning: Model has to be one of the following 'standard', 'null', 'homogeneous', 'heterogeneous', or 'selection' (the default)."
  if(!is.element(criterion, c("AIC", "BIC"))){
    criterion <- "AIC"
    print("Warning: selection criterion must be either 'AIC' or 'BIC'. Using the default AIC.")
  }
  if(missing(mc_S) || missing(Nf_S)){
    mod <- "standard"
    print("Warning: No mutant counts or final population size under stressful condition given. Using the standard model to infer the mutation rate under permissive conditions.")
  } 
  if(mod == "standard"){
    conv_input <- check_input(mc_UT, Nf_UT, plateff = plateff, fit_m = fit_m)
    if(conv_input[[1]]){
      res <- julia_call(
        "estimu",
        as.integer(conv_input[[2]]), conv_input[[3]], conv_input[[4]][1], conv_input[[5]][1],
        need_return = "R"
      )
      res <- list(res[[1]], res[[2]])
      print(paste0("Model used for inference: ", res[[2]]$model[1]))
    }
  }
  if(mod == "null"){
    conv_input_UT <- check_input(mc_UT, Nf_UT, plateff = plateff, fit_m = fit_m)
    conv_input_S <- check_input(mc_S, Nf_S)
    if(conv_input_UT[[1]]&&conv_input_S[[1]]){
      res <- julia_call(
        "estimu_0",
        conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], conv_input_UT[[5]],
        need_return = "R"
      )
      print(paste0("Model used for inference: ", res[[2]]$model[1]))
      res <- list(res[[1]], res[[2]])
    }
  }
  if(mod == "homogeneous"){
    conv_input_UT <- check_input(mc_UT, Nf_UT, plateff = plateff, fit_m = fit_m)
    conv_input_S <- check_input(mc_S, Nf_S)
    if(conv_input_UT[[1]]&&conv_input_S[[1]]){
      res <- julia_call(
        "estimu_hom",
        conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], conv_input_UT[[5]],
        need_return = "R"
      )
      print(paste0("Model used for inference: ", res[[2]]$model[1]))
      res <- list(res[[1]], res[[2]])
    }
  }
  if(mod == "heterogeneous"){
    conv_input_UT <- check_input(mc_UT, Nf_UT, plateff = plateff, fit_m = fit_m)
    conv_input_S <- check_input(mc_S, Nf_S, f_on = f_on, rel_div_on = rel_div_on)
    if(conv_input_UT[[1]]&&conv_input_S[[1]]){
      res <- julia_call(
        "estimu_het",
        conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], conv_input_S[[6]], conv_input_S[[7]], conv_input_UT[[5]],
        need_return = "R"
      )
      print(paste0("Model used for inference: ", res[[2]]$model[1]))
      res <- list(res[[1]], res[[2]])
    }
  }
  if(mod == "selection"){
    conv_input_UT <- check_input(mc_UT, Nf_UT, plateff = plateff, fit_m = fit_m)
    conv_input_S <- check_input(mc_S, Nf_S, rel_div_on = rel_div_on)
    if(conv_input_UT[[1]]&&conv_input_S[[1]]){
      if(conv_input_UT[[5]][1] != 1. && conv_input_UT[[5]][2] != 1.){
        res_null <- julia_call(
          "estimu_0",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], conv_input_UT[[5]],
          need_return = "R"
        )
        res_hom <- julia_call(
          "estimu_hom",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], conv_input_UT[[5]],
          need_return = "R"
        )
        hom <- list(res_null, res_hom)
        LRT_hom <- LRT(c(res_null[[2]]$LL[1], res_hom[[2]]$LL[1]))
        if (LRT_hom == 1){
          hom[[2]][[2]]$selection_result[1] <- "ns"
        } else {
          hom[[LRT_hom]][[2]]$selection_result[1] <- "best hom."
        }
        if(criterion == "AIC"){
          crit_hom <- hom[[LRT_hom]][[2]]$AIC[1] 
        } else {
          crit_hom <- hom[[LRT_hom]][[2]]$BIC[1] 
        }
        res_het_0 <- julia_call(
          "estimu_het",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, 0., conv_input_UT[[5]],
          need_return = "R"
        ) 
        if (conv_input_S[[7]] == 0.){
          res_het <- julia_call(
            "estimu_het",
            conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, FALSE, conv_input_UT[[5]],
            need_return = "R"
          ) 
          LRT_het <- min(LRT(c(res_null[[2]]$LL[1], res_het_0[[2]]$LL[1], -Inf, res_het[[2]]$LL[1])), 3)
        } else {
          res_het <- julia_call(
            "estimu_het",
            conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, conv_input_S[[7]], conv_input_UT[[5]],
            need_return = "R"
          ) 
          LRT_het <- LRT(c(res_null[[2]]$LL[1], res_het_0[[2]]$LL[1], res_het[[2]]$LL[1]))
        }
        het <- list(res_null, res_het_0, res_het)
        if (LRT_het == 1){
          het[[2]][[2]]$selection_result[1] <- "ns"
          het[[3]][[2]]$selection_result[1] <- "ns"
        } else {
          het[[LRT_het]][[2]]$selection_result[1] <- "best het."
        }
        if(criterion == "AIC"){
          crit_het <- het[[LRT_het]][[2]]$AIC[1] 
        } else {
          crit_het <- het[[LRT_het]][[2]]$BIC[1]
        }
        if(LRT_hom > 1 || LRT_het > 1){
          print("Significant support for SIM")
          if(crit_hom - crit_het < -2){
            hom[[LRT_hom]][[2]]$selection_result[1] <- "selected"
            print(paste0("Selected model: ", hom[[LRT_hom]][[2]]$model[1]))
          } else {
            if(crit_hom - crit_het > 2){
              het[[LRT_het]][[2]]$selection_result[1] <- "selected"
              print(paste0("Selected model: ", het[[LRT_het]][[2]]$model[1]))
            } else {
              print(paste0("Model selection between ", hom[[LRT_hom]][[2]]$model[1], "and", het[[LRT_het]][[2]]$model[1], " inconclusive."))
            }
          }
        } else {
          print("No significant support for SIM")
        }
        print("Estimated parameters under all models are:")
        res <- list(hom[[1]][[2]], hom[[1]][[1]], hom[[2]][[2]], hom[[2]][[1]], het[[2]][[2]], het[[2]][[1]], het[[3]][[2]], het[[3]][[1]])
      } else {
        res_null <- julia_call(
          "estimu_0",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]],
          need_return = "R"
        )
        res_null_constr <- julia_call(
          "estimu_0",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE,
          need_return = "R"
        )
        res_hom_1 <- julia_call(
          "estimu_hom",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]],
          need_return = "R"
        )
        res_hom_constr <- julia_call(
          "estimu_hom",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE,
          need_return = "R"
        )
        res_hom_unconstr <- julia_call(
          "estimu_hom",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], c(FALSE, FALSE),
          need_return = "R"
        )
        LRT_null <- LRT(c(res_null[[2]]$LL[1], res_null_constr[[2]]$LL[1]))
        hom <- list(res_null, res_hom_1, res_hom_constr, res_hom_unconstr)
        LRT_hom <- LRT(c(res_null[[2]]$LL[1], res_hom_1[[2]]$LL[1], res_hom_constr[[2]]$LL[1], res_hom_unconstr[[2]]$LL[1]))
        if (LRT_null == 2 && LRT_hom == 2){
          if(criterion == "AIC"){
            crit_null <- res_null_constr[[2]]$AIC[1]
            crit_hom <- res_hom_1[[2]]$AIC[1] 
          } else {
            crit_null <- res_null_constr[[2]]$BIC[1]
            crit_hom <- res_hom_1[[2]]$BIC[1] 
          }
          if(crit_null - crit_hom < -2){
            LRT_hom <- 1
          } else {
            if(crit_null - crit_hom < 2){
              print(paste0("Model selection between ", res_null_constr[[2]]$model[1], "and", res_hom_1[[2]]$model[1], " inconclusive."))
            }
          }
        }
        if (LRT_hom == 1){
          hom[[2]][[2]]$selection_result[1] <- "ns"
          hom[[3]][[2]]$selection_result[1] <- "ns"
          hom[[4]][[2]]$selection_result[1] <- "ns"
        } else {
          hom[[LRT_hom]][[2]]$selection_result[1] <- "best hom."
        }
        if(criterion == "AIC"){
          crit_hom <- hom[[LRT_hom]][[2]]$AIC[1] 
        } else {
          crit_hom <- hom[[LRT_hom]][[2]]$BIC[1] 
        }
        res_het_0 <- julia_call(
          "estimu_het",
          conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, 0.,
          need_return = "R"
        ) 
        if (conv_input_S[[7]] == 0.){
          res_het <- julia_call(
            "estimu_het",
            conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, FALSE,
            need_return = "R"
          ) 
          LRT_het <- min(LRT(c(res_null[[2]]$LL[1], res_het_0[[2]]$LL[1], -Inf, res_het[[2]]$LL[1])), 3)
        } else {
          res_het <- julia_call(
            "estimu_het",
            conv_input_UT[[2]], conv_input_UT[[3]], conv_input_S[[2]], conv_input_S[[3]], conv_input_UT[[4]], FALSE, conv_input_S[[7]],
            need_return = "R"
          )  
          LRT_het <- LRT(c(res_null[[2]]$LL[1], res_het_0[[2]]$LL[1], res_het[[2]]$LL[1]))
        }
        het <- list(res_null, res_het_0, res_het)
        if (LRT_het == 1){
          het[[2]][[2]]$selection_result[1] <- "ns"
          het[[3]][[2]]$selection_result[1] <- "ns"
        } else {
          het[[LRT_het]][[2]]$selection_result[1] <- "best het."
        }
        if(criterion == "AIC"){
          crit_het <- het[[LRT_het]][[2]]$AIC[1] 
        } else {
          crit_het <- het[[LRT_het]][[2]]$BIC[1]
        }
        if(LRT_hom > 1 || LRT_het > 1){
          print("Significant support for SIM")
          if(crit_hom - crit_het < -2){
            hom[[LRT_hom]][[2]]$selection_result[1] <- "selected"
            print(paste0("Selected model: ", hom[[LRT_hom]][[2]]$model[1]))
          } else {
            if(crit_hom - crit_het > 2){
              het[[LRT_het]][[2]]$selection_result[1] <- "selected"
              print(paste0("Selected model: ", het[[LRT_het]][[2]]$model[1]))
            } else {
              print(paste0("Model selection between ", hom[[LRT_hom]][[2]]$model[1], " and ", het[[LRT_het]][[2]]$model[1], " inconclusive."))
            }
          }
        } else {
          print("No significant support for SIM")
        }
        print("Estimated parameters under all models are:")
        res <- list(res_null[[2]], res_null[[1]], res_null_constr[[2]], res_null_constr[[1]], hom[[2]][[2]], hom[[2]][[1]], hom[[3]][[2]], hom[[3]][[1]], hom[[4]][[2]], hom[[4]][[1]], het[[2]][[2]], het[[2]][[1]], het[[3]][[2]], het[[3]][[1]])
      }
    }
  }
  return(res)
}

pMudi <- function(mc_max, Nf, mu, plateff=1, fit_m=1., S=FALSE, f_on=0., rel_div_on=0.){
  if(S == FALSE){
    p <- julia_call(
      "p_mudi",
      as.integer(mc_max), Nf, mu, fit_m, plateff,
      need_return = "R"
    )
  } else {
    p <- julia_call(
      "p_mudi",
      as.integer(mc_max), Nf, mu, S, f_on, rel_div_on, fit_m, plateff,
      need_return = "R"
    )
  }
  return(p)
}
  
read_counts <- function(df_row){
  v1 <- as.vector(df_row)
  v2 <- c()
  for (j in 1:length(v1)){
    if(!is.na(as.numeric(v1[[j]]))){
      v2 <- append(v2, as.numeric(v1[[j]]))
    }
  }
  return(v2)
}
