meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]
est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]

chisq_1_95 <- 3.84145882069412447634704221854917705059051513671875
chisq_2_95 <- 5.99146454710797993215010137646459043025970458984375
chisq_3_95 <- 7.81472790325117738774451936478726565837860107421875

hom_by_LRT <- c()
M <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_UT <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_ratio <- matrix(nrow = length(meta_data$ID), ncol = 3)
hom_by_LRT_constr <- c()
M_constr <- matrix(nrow = length(meta_data$ID), ncol = 3)
fitm_UT_constr <- matrix(nrow = length(meta_data$ID), ncol = 3)
het_by_LRT <- c()
S <- matrix(nrow = length(meta_data$ID), ncol = 3)
rel_div_on <- matrix(nrow = length(meta_data$ID), ncol = 3)
SIM <- c()
SIM_constr <- c()
by_AIC <- c()
by_AIC_constr <- c()
by_BIC <- c()
by_BIC_constr <- c()

for (i in 1:length(meta_data$ID)) {
  df_i <- subset(subset(est_paras, is.element(model, c("no_SIM", "hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr"))), ID == meta_data$ID[i])
  s <- 1
  s_constr <- 1
  if (df_i$LL[1] + chisq_1_95/2 < df_i$LL[2]) {
    s <- 2
    s_constr <- 2
    if (df_i$LL[2] + chisq_1_95/2 < df_i$LL[3]) {
      s <- 3
      s_constr <- 3
      if (df_i$LL[3] + chisq_1_95/2 < df_i$LL[4]) {
        s <- 4
      }        
    } else {
      if (df_i$LL[2] + chisq_2_95/2 < df_i$LL[4]){
        s <- 4
      }
    }    
  } else {
    if (df_i$LL[1] + chisq_2_95/2 < df_i$LL[3]) {
      s <- 3
      s_constr <- 3
      if (df_i$LL[3] +chisq_1_95/2 < df_i$LL[4]) {
        s <- 4
      }
    } else {
      if (df_i$LL[1] + chisq_3_95/2 < df_i$LL[4]){
        s <- 4
      }
    }
  }
  M[i,1] <- df_i$M_MLE[s]
  M[i,2] <- df_i$M_lower_bound[s]
  M[i,3] <- df_i$M_upper_bound[s]
  fitm_UT[i,1] <- df_i$fitm_UT_MLE[s]
  fitm_UT[i,2] <- df_i$fitm_UT_lower_bound[s]
  fitm_UT[i,3] <- df_i$fitm_UT_upper_bound[s]
  fitm_ratio[i,1] <- df_i$fitm_ratio_MLE[s]
  fitm_ratio[i,2] <- df_i$fitm_ratio_lower_bound[s]
  fitm_ratio[i,3] <- df_i$fitm_ratio_upper_bound[s]
  hom_by_LRT[i] <- c("no_SIM", "hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr")[s]
  AIC_hom <- df_i$AIC[s]
  BIC_hom <- df_i$BIC[s]
  M_constr[i,1] <- df_i$M_MLE[s_constr]
  M_constr[i,2] <- df_i$M_lower_bound[s_constr]
  M_constr[i,3] <- df_i$M_upper_bound[s_constr]
  fitm_UT_constr[i,1] <- df_i$fitm_UT_MLE[s]
  fitm_UT_constr[i,2] <- df_i$fitm_UT_lower_bound[s]
  fitm_UT_constr[i,3] <- df_i$fitm_UT_upper_bound[s]
  hom_by_LRT_constr[i] <- c("no_SIM", "hom_wo_fitm", "hom_fitm")[s_constr]
  AIC_hom_constr <- df_i$AIC[s_constr]
  BIC_hom_constr <- df_i$BIC[s_constr]
  df_i <- subset(subset(est_paras, is.element(model, c("no_SIM", "het_zero_div", "het_div"))), ID == meta_data$ID[i])
  s <- 1
  if (df_i$LL[1] + chisq_1_95/2 < df_i$LL[2]) {
    s <- 2
    if (df_i$LL[2] + chisq_1_95/2 < df_i$LL[3]) {
      s <- 3
    }
  } else {
    if (df_i$LL[1] + chisq_3_95/2 < df_i$LL[3]) {
      s <- 3
    }
  }
  S[i,1] <- df_i$S_MLE[s]
  S[i,2] <- df_i$S_lower_bound[s]
  S[i,3] <- df_i$S_upper_bound[s]
  rel_div_on[i,1] <- df_i$rel_div_on_MLE[s]
  rel_div_on[i,2] <- df_i$rel_div_on_lower_bound[s]
  rel_div_on[i,3] <- df_i$rel_div_on_upper_bound[s]
  het_by_LRT[i] <- c("no_SIM", "het_zero_div", "het_div")[s]
  AIC_het <- df_i$AIC[s]
  BIC_het <- df_i$BIC[s]
  if (het_by_LRT[i] == "no_SIM") {
    if (hom_by_LRT[i] == "no_SIM") {
      SIM[i] <- "no"
      by_AIC[i] <- "-"
      by_BIC[i] <- "-"
    } else {
      SIM[i] <- "hom" 
      if (M[i,1] <= 1) {
        SIM[i] <- "no" 
      }
      by_AIC[i] <- "-"
      by_BIC[i] <- "-"
    }
    if (hom_by_LRT_constr[i] == "no_SIM") {
      SIM_constr[i] <- "no"
      by_AIC_constr[i] <- "-"
      by_BIC_constr[i] <- "-"
    } else {
      SIM_constr[i] <- "hom"
      if (M_constr[i,1] <= 1) {
        SIM_constr[i] <- "no" 
      }
      by_AIC_constr[i] <- "-"
      by_BIC_constr[i] <- "-"
    }
  } else {
    if (hom_by_LRT[i] == "no_SIM") {
      SIM[i] <- "het"
      by_AIC[i] <- "-"
      by_BIC[i] <- "-"
    } else {
      SIM[i] <- "no"
      by_AIC[i] <- "-"
      by_BIC[i] <- "-"
      if (M[i,1] > 1){
        SIM[i] <- "yes"
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
      }
    }
    if (hom_by_LRT_constr[i] == "no_SIM") {
      SIM_constr[i] <- "het"
      by_AIC_constr[i] <- "-"
      by_BIC_constr[i] <- "-"
    } else {
      SIM_constr[i] <- "no"
      by_AIC_constr[i] <- "-"
      by_BIC_constr[i] <- "-"
      if (M_constr[i,1] > 1){
        SIM_constr[i] <- "yes"
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
}

est_sum <- data.frame(ID=meta_data$ID, hom_by_LRT=hom_by_LRT, M=M, fitm_UT=fitm_UT, fitm_ratio=fitm_ratio, hom_by_LRT_constr=hom_by_LRT_constr, M_constr=M_constr, fitm_UT_constr=fitm_UT_constr, het_by_LRT=het_by_LRT, S=S, rel_div_on=rel_div_on, SIM=SIM, SIM_constr=SIM_constr, by_AIC=by_AIC, by_AIC_constr=by_AIC_constr, by_BIC=by_BIC, by_BIC_constr=by_BIC_constr)
write.csv(est_sum, file = "experimental_data/est_sum.csv")
