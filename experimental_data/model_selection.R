library("plyr")

est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
IDs <- unique(est_paras$ID)
est_sum <- data.frame(ID = IDs)

est_sum$p_value_UT_max <- pmax(
  subset(est_paras, model=="hom_wo_fitm")$p_value_UT, 
  subset(est_paras, model=="hom_fitm")$p_value_UT, 
  subset(est_paras, model=="hom_fitm_unconstr")$p_value_UT)

est_sum$p_value_test_min <- pmin(
  subset(est_paras, model=="hom_wo_fitm")$p_value_test, 
  subset(est_paras, model=="hom_fitm_unconstr")$p_value_test,
  na.rm = TRUE)

est_sum$M_wo_fitm <- cbind(
  subset(est_paras, model=="hom_wo_fitm")$M_MLE, 
  subset(est_paras, model=="hom_wo_fitm")$M_lower_bound, 
  subset(est_paras, model=="hom_wo_fitm")$M_upper_bound)

m <- c("no_SIM_wo_fitm","no_SIM_fitm","no_SIM_fitm_unconstr","hom_wo_fitm","hom_fitm","hom_fitm_unconstr")
est_paras_hom_only <- subset(est_paras, is.element(model, m))

min_AIC_hom <- pmin(
  subset(est_paras_hom_only, model==m[1])$AIC_joint, 
  subset(est_paras_hom_only, model==m[2])$AIC_joint, 
  subset(est_paras_hom_only, model==m[3])$AIC_joint, 
  subset(est_paras_hom_only, model==m[4])$AIC_joint, 
  subset(est_paras_hom_only, model==m[5])$AIC_joint, 
  subset(est_paras_hom_only, model==m[6])$AIC_joint)
est_paras_hom_AIC <- est_paras_hom_only[est_paras_hom_only$AIC_joint == rep(min_AIC_hom, each = length(m)), ]
est_sum$hom_by_AIC <- est_paras_hom_AIC$model
est_sum$M_AIC <- cbind(
  est_paras_hom_AIC$M_MLE,
  est_paras_hom_AIC$M_lower_bound,
  est_paras_hom_AIC$M_upper_bound)
est_sum$fitm_ratio_AIC <- cbind(
  est_paras_hom_AIC$fitm_ratio_MLE,
  est_paras_hom_AIC$fitm_ratio_lower_bound,
  est_paras_hom_AIC$fitm_ratio_upper_bound)

min_BIC_hom <- pmin(
  subset(est_paras_hom_only, model==m[1])$BIC_joint, 
  subset(est_paras_hom_only, model==m[2])$BIC_joint, 
  subset(est_paras_hom_only, model==m[3])$BIC_joint, 
  subset(est_paras_hom_only, model==m[4])$BIC_joint, 
  subset(est_paras_hom_only, model==m[5])$BIC_joint, 
  subset(est_paras_hom_only, model==m[6])$BIC_joint)
est_paras_hom_BIC <- est_paras_hom_only[est_paras_hom_only$BIC_joint == rep(min_BIC_hom, each = length(m)), ]
est_sum$hom_by_BIC <- est_paras_hom_BIC$model
est_sum$M_BIC <- cbind(
  est_paras_hom_BIC$M_MLE,
  est_paras_hom_BIC$M_lower_bound,
  est_paras_hom_BIC$M_upper_bound)
est_sum$fitm_ratio_BIC <- cbind(
  est_paras_hom_BIC$fitm_ratio_MLE,
  est_paras_hom_BIC$fitm_ratio_lower_bound,
  est_paras_hom_BIC$fitm_ratio_upper_bound)

m <- c("het_zero_div","het_div","het_div_fon")
est_paras_het_only <- subset(est_paras, is.element(model, m))

min_AIC_het <- pmin(
  subset(est_paras_het_only, model==m[1])$AIC_joint, 
  subset(est_paras_het_only, model==m[2])$AIC_joint, 
  subset(est_paras_het_only, model==m[3])$AIC_joint)
est_paras_het_AIC <- est_paras_het_only[est_paras_het_only$AIC_joint == rep(min_AIC_het, each = length(m)), ]
est_sum$het_by_AIC <- est_paras_het_AIC$model
est_sum$S_AIC <- cbind(
  est_paras_het_AIC$S_MLE,
  est_paras_het_AIC$S_lower_bound,
  est_paras_het_AIC$S_upper_bound)
est_sum$rel_div_on_AIC <- cbind(
  est_paras_het_AIC$rel_div_on_MLE,
  est_paras_het_AIC$rel_div_on_lower_bound,
  est_paras_het_AIC$rel_div_on_upper_bound)

min_BIC_het <- pmin(
  subset(est_paras_het_only, model==m[1])$BIC_joint, 
  subset(est_paras_het_only, model==m[2])$BIC_joint, 
  subset(est_paras_het_only, model==m[3])$BIC_joint)
est_paras_het_BIC <- est_paras_het_only[est_paras_het_only$BIC_joint == rep(min_BIC_het, each = length(m)), ]
est_sum$het_by_BIC <- est_paras_het_BIC$model
est_sum$rel_div_on_BIC <- cbind(
  est_paras_het_BIC$rel_div_on_MLE,
  est_paras_het_BIC$rel_div_on_lower_bound,
  est_paras_het_BIC$rel_div_on_upper_bound)

m <- c("no_SIM_wo_fitm","no_SIM_fitm","no_SIM_fitm_unconstr","hom_wo_fitm","hom_fitm","hom_fitm_unconstr","het_zero_div","het_div","het_div_fon")
est_paras_all <- subset(est_paras, is.element(model, m))

min_AIC <- pmin(min_AIC_hom, min_AIC_het)
est_paras_AIC <- est_paras_all[est_paras_all$AIC_joint == rep(min_AIC, each = length(m)), ]
est_sum$by_AIC <- est_paras_AIC$model

min_BIC <- pmin(min_BIC_hom, min_BIC_het)
est_paras_BIC <- est_paras_all[est_paras_all$BIC_joint == rep(min_BIC, each = length(m)), ]
est_sum$by_BIC <- est_paras_BIC$model

write.csv(est_sum, file = "experimental_data/est_sum.csv")
