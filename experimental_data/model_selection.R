library("plyr")

meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]
est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]

LRT <- function(LLs){
  q <- numeric()
  for (d in 1:(length(LLs)-1)) {
    q <- append(q, qchisq(0.95, df=d)/2)
  }  
  m <- 1
  while (m < length(LLs)) {
    for (d in 1:(length(LLs)-m)) {
      if (LLs[m] + q[d] < LLs[m+d]) {
        m <- m + d
        break
      }      
    }
    if (d == (length(LLs)-m)) {
      break
    }
  }
  return(m)
}

M_wo_fitm <- matrix(nrow = length(meta_data$ID), ncol = 3)
SIM <- logical(length(meta_data$ID))
hom_by_LRT <- character(length(meta_data$ID))
M <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_UT <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_ratio <- matrix(nrow = length(meta_data$ID), ncol = 3)
hom_by_LRT_constr <- character(length(meta_data$ID))
M_constr <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_UT_constr <- matrix(nrow = length(meta_data$ID), ncol = 3)
het_by_LRT <- character(length(meta_data$ID))
S <- matrix(nrow = length(meta_data$ID), ncol = 3)
rel_div_on <- matrix(nrow = length(meta_data$ID), ncol = 3)
M_het <- matrix(nrow = length(meta_data$ID), ncol = 3)
by_AIC <- character(length(meta_data$ID))
by_AIC_constr <- character(length(meta_data$ID))
by_BIC <- character(length(meta_data$ID))
by_BIC_constr <- character(length(meta_data$ID))

for (i in 1:length(meta_data$ID)) {
  df_SIM <- subset(subset(est_paras, is.element(model, c("no_SIM", "hom_wo_fitm"))), ID == meta_data$ID[i])
  M_wo_fitm[i,1] <- df_SIM$M_MLE[2]
  M_wo_fitm[i,2] <- df_SIM$M_lower_bound[2]
  M_wo_fitm[i,3] <- df_SIM$M_upper_bound[2]
  if ((LRT(df_SIM$LL) == 2) && M_wo_fitm[i,1] > 1) {
    SIM[i] <- TRUE
    df_hom <- subset(subset(est_paras, is.element(model, c("hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr"))), ID == meta_data$ID[i])
    s <- LRT(df_hom$LL)
    s_constr <- LRT(df_hom$LL[1:2])
    M[i,1] <- df_hom$M_MLE[s]
    M[i,2] <- df_hom$M_lower_bound[s]
    M[i,3] <- df_hom$M_upper_bound[s]
    fitm_UT[i,1] <- df_hom$fitm_UT_MLE[s]
    fitm_UT[i,2] <- df_hom$fitm_UT_lower_bound[s]
    fitm_UT[i,3] <- df_hom$fitm_UT_upper_bound[s]
    fitm_ratio[i,1] <- df_hom$fitm_ratio_MLE[s]
    fitm_ratio[i,2] <- df_hom$fitm_ratio_lower_bound[s]
    fitm_ratio[i,3] <- df_hom$fitm_ratio_upper_bound[s]
    hom_by_LRT[i] <- c("hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr")[s]
    AIC_hom <- df_hom$AIC[s]
    BIC_hom <- df_hom$BIC[s]
    M_constr[i,1] <- df_hom$M_MLE[s_constr]
    M_constr[i,2] <- df_hom$M_lower_bound[s_constr]
    M_constr[i,3] <- df_hom$M_upper_bound[s_constr]
    fitm_UT_constr[i,1] <- df_hom$fitm_UT_MLE[s]
    fitm_UT_constr[i,2] <- df_hom$fitm_UT_lower_bound[s]
    fitm_UT_constr[i,3] <- df_hom$fitm_UT_upper_bound[s]
    hom_by_LRT_constr[i] <- c("hom_wo_fitm", "hom_fitm")[s_constr]
    AIC_hom_constr <- df_hom$AIC[s_constr]
    BIC_hom_constr <- df_hom$BIC[s_constr] 
    df_het <- subset(subset(est_paras, is.element(model, c("no_SIM", "het_zero_div", "het_div_fon", "het_div"))), ID == meta_data$ID[i])
    df_het <- arrange(df_het, match(model, c("no_SIM", "het_zero_div", "het_div_fon", "het_div")))
    s <- LRT(df_het$LL)
    if (s == 1) {
      by_AIC[i] <- "hom"
      by_BIC[i] <- "hom"
      by_AIC_constr[i] <- "hom"
      by_BIC_constr[i] <- "hom"
    } else {
      S[i,1] <- df_het$S_MLE[s]
      S[i,2] <- df_het$S_lower_bound[s]
      S[i,3] <- df_het$S_upper_bound[s]
      rel_div_on[i,1] <- df_het$rel_div_on_MLE[s]
      rel_div_on[i,2] <- df_het$rel_div_on_lower_bound[s]
      rel_div_on[i,3] <- df_het$rel_div_on_upper_bound[s]
      M_het[i,1] <- df_het$M_MLE[s]
      M_het[i,2] <- df_het$M_lower_bound[s]
      M_het[i,3] <- df_het$M_upper_bound[s]
      het_by_LRT[i] <- c("no_SIM", "het_zero_div", "het_div_fon", "het_div")[s]
      AIC_het <- df_het$AIC[s]
      BIC_het <- df_het$BIC[s] 
      if (AIC_het - AIC_hom < -2){
        by_AIC[i] <- "het"
      } else {
        if (AIC_het - AIC_hom >2){
          by_AIC[i] <- "hom"
        } else {
          by_AIC[i] <- "none"
        }
      }
      if (BIC_het - BIC_hom < -2){
        by_BIC[i] <- "het"
      } else {
        if (BIC_het - BIC_hom >2){
          by_BIC[i] <- "hom"
        } else {
          by_BIC[i] <- "none"
        }
      }
      if (AIC_het - AIC_hom_constr > 2){
        by_AIC_constr[i] <- "hom"
      } else {
        if (AIC_het - AIC_hom_constr < -2){
          by_AIC_constr[i] <- "het"
        } else {
          by_AIC_constr[i] <- "none"
        }
      }
      if (BIC_het - BIC_hom_constr > 2){
        by_BIC_constr[i] <- "hom"
      } else {
        if (BIC_het - BIC_hom_constr < -2){
          by_BIC_constr[i] <- "het"
        } else {
          by_BIC_constr[i] <- "none"
        }
      }
    }
  }
}

est_sum <- data.frame(ID=meta_data$ID, M_wo_fitm=M_wo_fitm, SIM=SIM, hom_by_LRT=hom_by_LRT, M=M, fitm_UT=fitm_UT, fitm_ratio=fitm_ratio, hom_by_LRT_constr=hom_by_LRT_constr, M_constr=M_constr, fitm_UT_constr=fitm_UT_constr, het_by_LRT=het_by_LRT, S=S, rel_div_on=rel_div_on, M_het=M_het, by_AIC=by_AIC, by_AIC_constr=by_AIC_constr, by_BIC=by_BIC, by_BIC_constr=by_BIC_constr)
write.csv(est_sum, file = "experimental_data/est_sum.csv")
