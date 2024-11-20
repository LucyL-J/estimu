library("plyr")

est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
IDs <- unique(est_paras$ID)

M_wo_fitm <- matrix(nrow = length(IDs), ncol = 3)
SIM <- logical(length(IDs))
hom_by_LRT <- character(length(IDs))
M <- matrix(nrow = length(IDs), ncol = 3)
fitm_UT <- matrix(nrow = length(IDs), ncol = 3)
fitm_ratio <- matrix(nrow = length(IDs), ncol = 3)
hom_by_LRT_constr <- character(length(IDs))
M_constr <- matrix(nrow = length(IDs), ncol = 3)
fitm_UT_constr <- matrix(nrow = length(IDs), ncol = 3)
het_by_LRT <- character(length(IDs))
S <- matrix(nrow = length(IDs), ncol = 3)
rel_div_on <- matrix(nrow = length(IDs), ncol = 3)
M_het <- matrix(nrow = length(IDs), ncol = 3)
by_AIC <- character(length(IDs))
Delta_AIC <- numeric(length(IDs))
by_AIC_constr <- character(length(IDs))
Delta_AIC_constr <- numeric(length(IDs))
by_BIC <- character(length(IDs))
Delta_BIC <- numeric(length(IDs))
by_BIC_constr <- character(length(IDs))
Delta_BIC_constr <- numeric(length(IDs))

for (i in 1:length(IDs)) {
  df_SIM <- subset(subset(est_paras, is.element(model, c("no_SIM_wo_fitm", "hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr"))), ID == IDs[i])
  M_wo_fitm[i,1] <- df_SIM$M_MLE[2]
  M_wo_fitm[i,2] <- df_SIM$M_lower_bound[2]
  M_wo_fitm[i,3] <- df_SIM$M_upper_bound[2]
  s <- LRT(df_SIM$LL)
  df_null <- subset(subset(est_paras, is.element(model, c("no_SIM_wo_fitm", "no_SIM_fitm"))), ID == IDs[i])
  s_null <- LRT(df_null$LL)
  if ((s >= s_null) && (df_SIM$AIC[s] - df_null$AIC[s_null] < -4) && (M_wo_fitm[i,1] > 1)) {
    SIM[i] <- TRUE
    df_hom <- subset(subset(est_paras, is.element(model, c("hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr"))), ID == IDs[i])
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
    fitm_UT_constr[i,1] <- df_hom$fitm_UT_MLE[s_constr]
    fitm_UT_constr[i,2] <- df_hom$fitm_UT_lower_bound[s_constr]
    fitm_UT_constr[i,3] <- df_hom$fitm_UT_upper_bound[s_constr]
    hom_by_LRT_constr[i] <- c("hom_wo_fitm", "hom_fitm")[s_constr]
    AIC_hom_constr <- df_hom$AIC[s_constr]
    BIC_hom_constr <- df_hom$BIC[s_constr] 
    df_het <- subset(subset(est_paras, is.element(model, c("no_SIM_wo_fitm", "het_zero_div", "het_div_fon", "het_div"))), ID == IDs[i])
    df_het <- arrange(df_het, match(model, c("no_SIM_wo_fitm", "het_zero_div", "het_div_fon", "het_div")))
    s <- LRT(df_het$LL)
    if (s == 1) {
      by_AIC[i] <- "hom"
      by_BIC[i] <- "hom"
      by_AIC_constr[i] <- "hom"
      by_BIC_constr[i] <- "hom"
      Delta_AIC[i] <- Inf
      Delta_BIC[i] <- Inf
      Delta_AIC_constr[i] <- Inf
      Delta_BIC_constr[i] <- Inf
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
      het_by_LRT[i] <- c("no_SIM_wo_fitm", "het_zero_div", "het_div_fon", "het_div")[s]
      AIC_het <- df_het$AIC[s]
      BIC_het <- df_het$BIC[s] 
      Delta_AIC[i] <- AIC_het - AIC_hom
      Delta_BIC[i] <- BIC_het - BIC_hom
      Delta_AIC_constr[i] <- AIC_het - AIC_hom_constr
      Delta_BIC_constr[i] <- BIC_het - BIC_hom_constr
      if (AIC_het - AIC_hom < -4){
        by_AIC[i] <- "het"
      } else {
        if (AIC_het - AIC_hom > 4){
          by_AIC[i] <- "hom"
        } else {
          by_AIC[i] <- "none"
        }
      }
      if (BIC_het - BIC_hom < -4){
        by_BIC[i] <- "het"
      } else {
        if (BIC_het - BIC_hom > 4){
          by_BIC[i] <- "hom"
        } else {
          by_BIC[i] <- "none"
        }
      }
      if (AIC_het - AIC_hom_constr > 4){
        by_AIC_constr[i] <- "hom"
      } else {
        if (AIC_het - AIC_hom_constr < -4){
          by_AIC_constr[i] <- "het"
        } else {
          by_AIC_constr[i] <- "none"
        }
      }
      if (BIC_het - BIC_hom_constr > 4){
        by_BIC_constr[i] <- "hom"
      } else {
        if (BIC_het - BIC_hom_constr < -4){
          by_BIC_constr[i] <- "het"
        } else {
          by_BIC_constr[i] <- "none"
        }
      }
    }
  }
}

est_sum <- data.frame(ID=IDs, M_wo_fitm=M_wo_fitm, SIM=SIM, hom_by_LRT=hom_by_LRT, M=M, fitm_UT=fitm_UT, fitm_ratio=fitm_ratio, hom_by_LRT_constr=hom_by_LRT_constr, M_constr=M_constr, fitm_UT_constr=fitm_UT_constr, het_by_LRT=het_by_LRT, S=S, rel_div_on=rel_div_on, M_het=M_het, by_AIC=by_AIC, Delta_AIC=Delta_AIC, by_AIC_constr=by_AIC_constr, Delta_AIC_constr=Delta_AIC_constr, by_BIC=by_BIC, Delta_BIC=Delta_BIC, by_BIC_constr=by_BIC_constr, Delta_BIC_constr=Delta_BIC_constr)
write.csv(est_sum, file = "experimental_data/est_sum.csv")
