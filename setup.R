check_input <- function(mc, Nf, eff=1, fit_m=1, rel_div_on=FALSE, f_on=0.1){
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
  if(is.numeric(eff) && (sum(eff>0)==length(eff)) && (sum(eff<=1)==length(eff))){
    if(length(eff)>2){
      eff <- eff[1:2]
      print("Note: More than two values for plating efficency given; only the first two (untreated & stressful) will be used in the inference.")
    }
  } else {
    status <- FALSE
    print("Error: Plating efficency(s) must be strictly between zero and one.")
  }
  if(is.numeric(fit_m) && (sum(fit_m>=0)==length(fit_m))){
    if(length(fit_m)>2){
      fit_m <- fit_m[1:2]
      print("Note: More than two values for mutant fitness given; only the first two (untreated & stressful) will be used in the inference.")
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
  return(list(status, mc, as.numeric(Nf), fit_m, f_on, rel_div_on))
}
