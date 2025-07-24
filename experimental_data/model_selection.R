library("plyr")

est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
IDs <- unique(est_paras$ID)
est_sum <- data.frame(ID = IDs)

est_sum$p_value_UT_max <- pmax(
  subset(est_paras, model=="HOM0")$p_value_UT, 
  subset(est_paras, model=="HOM1")$p_value_UT, 
  subset(est_paras, model=="HOM2")$p_value_UT)

est_sum$p_value_test_min <- pmin(
  subset(est_paras, model=="HOM0")$p_value_test, 
  subset(est_paras, model=="HOM2")$p_value_test,
  na.rm = TRUE)

est_sum$M_HOM0 <- cbind(
  subset(est_paras, model=="HOM0")$M_MLE, 
  subset(est_paras, model=="HOM0")$M_lower_bound, 
  subset(est_paras, model=="HOM0")$M_upper_bound)

m_null <- c("N0", "N1", "N2")
est_paras_null <- subset(est_paras, is.element(model, m_null))

min_AIC_null <- pmin(
  subset(est_paras_null, model==m_null[1])$AIC, 
  subset(est_paras_null, model==m_null[2])$AIC, 
  subset(est_paras_null, model==m_null[3])$AIC)
est_paras_null_AIC <- est_paras_null[est_paras_null$AIC == rep(min_AIC_null, each = length(m_null)), ]
est_sum$null_by_AIC <- est_paras_null_AIC$model

min_AIC_corr_null <- pmin(
  subset(est_paras_null, model==m_null[1])$AIC_corr, 
  subset(est_paras_null, model==m_null[2])$AIC_corr, 
  subset(est_paras_null, model==m_null[3])$AIC_corr)
est_paras_null_AIC_corr <- est_paras_null[est_paras_null$AIC_corr == rep(min_AIC_corr_null, each = length(m_null)), ]
est_sum$null_by_AIC_corr <- est_paras_null_AIC_corr$model

m_hom <- c("HOM0", "HOM1", "HOM2")
est_paras_hom <- subset(est_paras, is.element(model, m_hom))

min_AIC_hom <- pmin(
  subset(est_paras_hom, model==m_hom[1])$AIC, 
  subset(est_paras_hom, model==m_hom[2])$AIC, 
  subset(est_paras_hom, model==m_hom[3])$AIC)
est_paras_hom_AIC <- est_paras_hom[est_paras_hom$AIC == rep(min_AIC_hom, each = length(m_hom)), ]
est_sum$null_by_AIC <- est_paras_hom_AIC$model

min_AIC_corr_hom <- pmin(
  subset(est_paras_hom, model==m_hom[1])$AIC_corr, 
  subset(est_paras_hom, model==m_hom[2])$AIC_corr, 
  subset(est_paras_hom, model==m_hom[3])$AIC_corr)
est_paras_hom_AIC_corr <- est_paras_hom[est_paras_hom$AIC_corr == rep(min_AIC_corr_hom, each = length(m_hom)), ]
est_sum$null_by_AIC_corr <- est_paras_hom_AIC_corr$model

m_null_hom <- c("N0","N1","N2","HOM0","HOM1","HOM2")
est_paras_null_hom <- subset(est_paras, is.element(model, m_null_hom))

min_AIC_NH <- pmin(min_AIC_null, min_AIC_hom)
est_paras_AIC_NH <- est_paras_null_hom[est_paras_null_hom$AIC == rep(min_AIC_NH, each = length(m_null_hom)), ]
est_sum$by_AIC_NH <- est_paras_AIC_NH$model
est_sum$Delta_AIC_NH <- min_AIC_null - min_AIC_null_hom

min_AIC_corr_NH <- pmin(min_AIC_corr_null, min_AIC_corr_hom)
est_paras_AIC_corr_NH <- est_paras_null_hom[est_paras_null_hom$AIC_corr == rep(min_AIC_corr_NH, each = length(m_null_hom)), ]
est_sum$by_AIC_corr_NH <- est_paras_AIC_corr_NH$model
est_sum$Delta_AIC_corr_NH <- min_AIC_corr_null - min_AIC_corr_hom

min_AIC_null_hom <- pmin(
  subset(est_paras_null_hom, model==m_null_hom[1])$AIC, 
  subset(est_paras_null_hom, model==m_null_hom[2])$AIC, 
  subset(est_paras_null_hom, model==m_null_hom[3])$AIC, 
  subset(est_paras_null_hom, model==m_null_hom[4])$AIC, 
  subset(est_paras_null_hom, model==m_null_hom[5])$AIC, 
  subset(est_paras_null_hom, model==m_null_hom[6])$AIC)
est_paras_null_hom_AIC <- est_paras_null_hom[est_paras_null_hom$AIC == rep(min_AIC_null_hom, each = length(m_null_hom)), ]
est_sum$null_hom_by_AIC <- est_paras_null_hom_AIC$model
est_sum$M_AIC <- cbind(
  est_paras_null_hom_AIC$M_MLE,
  est_paras_null_hom_AIC$M_lower_bound,
  est_paras_null_hom_AIC$M_upper_bound)
est_sum$p_value_UT_hom_AIC <- est_paras_null_hom_AIC$p_value_UT
est_sum$p_value_S_hom_AIC <- est_paras_null_hom_AIC$p_value_S

min_AIC_corr_null_hom <- pmin(
  subset(est_paras_null_hom, model==m_null_hom[1])$AIC_corr, 
  subset(est_paras_null_hom, model==m_null_hom[2])$AIC_corr, 
  subset(est_paras_null_hom, model==m_null_hom[3])$AIC_corr, 
  subset(est_paras_null_hom, model==m_null_hom[4])$AIC_corr, 
  subset(est_paras_null_hom, model==m_null_hom[5])$AIC_corr, 
  subset(est_paras_null_hom, model==m_null_hom[6])$AIC_corr)
est_paras_null_hom_AIC_corr <- est_paras_null_hom[est_paras_null_hom$AIC_corr == rep(min_AIC_corr_null_hom, each = length(m_null_hom)), ]
est_sum$null_hom_by_AIC_corr <- est_paras_null_hom_AIC_corr$model
est_sum$M_AIC_corr <- cbind(
  est_paras_null_hom_AIC_corr$M_MLE,
  est_paras_null_hom_AIC_corr$M_lower_bound,
  est_paras_null_hom_AIC_corr$M_upper_bound)
est_sum$p_value_UT_hom_AIC_corr <- est_paras_null_hom_AIC_corr$p_value_UT
est_sum$p_value_S_hom_AIC_corr <- est_paras_null_hom_AIC_corr$p_value_S

m_null_het <- c("N0","HET0","HET2")
est_paras_het_only <- subset(est_paras, is.element(model, m_null_het))

min_AIC_null_het <- pmin(
  subset(est_paras_het_only, model==m_null_het[1])$AIC,
  subset(est_paras_het_only, model==m_null_het[2])$AIC,
  subset(est_paras_het_only, model==m_null_het[3])$AIC)
est_paras_null_het_AIC <- est_paras_het_only[est_paras_het_only$AIC == rep(min_AIC_null_het, each = length(m_null_het)), ]
est_sum$null_het_by_AIC <- est_paras_null_het_AIC$model
est_sum$S_AIC <- cbind(
  est_paras_null_het_AIC$S_MLE,
  est_paras_null_het_AIC$S_lower_bound,
  est_paras_null_het_AIC$S_upper_bound)
est_sum$p_value_UT_het_AIC <- est_paras_null_het_AIC$p_value_UT
est_sum$p_value_S_het_AIC <- est_paras_null_het_AIC$p_value_S

min_AIC_corr_null_het <- pmin(
  subset(est_paras_het_only, model==m_null_het[1])$AIC_corr, 
  subset(est_paras_het_only, model==m_null_het[2])$AIC_corr,
  subset(est_paras_het_only, model==m_null_het[3])$AIC_corr)
est_paras_null_het_AIC_corr <- est_paras_het_only[est_paras_het_only$AIC_corr == rep(min_AIC_corr_null_het, each = length(m_null_het)), ]
est_sum$null_het_by_AIC_corr <- est_paras_null_het_AIC_corr$model
est_sum$S_AIC_corr <- cbind(
  est_paras_null_het_AIC_corr$S_MLE,
  est_paras_null_het_AIC_corr$S_lower_bound,
  est_paras_null_het_AIC_corr$S_upper_bound)
est_sum$p_value_UT_het_AIC_corr <- est_paras_null_het_AIC_corr$p_value_UT
est_sum$p_value_S_het_AIC_corr <- est_paras_null_het_AIC_corr$p_value_S

m_all <- c("N0","N1","N2","HOM0","HOM1","HOM2","HET0","HET2")
est_paras_all <- subset(est_paras, is.element(model, m_all))

min_AIC <- pmin(min_AIC_null_hom, min_AIC_null_het)
est_paras_AIC <- est_paras_all[est_paras_all$AIC == rep(min_AIC, each = length(m_all)), ]
est_sum$by_AIC <- est_paras_AIC$model
est_sum$Delta_AIC <- min_AIC_null_hom - min_AIC_null_het

min_AIC_corr <- pmin(min_AIC_corr_null_hom, min_AIC_corr_null_het)
est_paras_AIC_corr <- est_paras_all[est_paras_all$AIC_corr == rep(min_AIC_corr, each = length(m_all)), ]
est_sum$by_AIC_corr <- est_paras_AIC_corr$model
est_sum$Delta_AIC_corr <- min_AIC_corr_null_hom - min_AIC_corr_null_het

write.csv(est_sum, file = "experimental_data/est_sum.csv")

for (id in IDs) {
  msel_res <- data.frame(model = m_all)
  est_paras_ID <- subset(est_paras_all, ID == id)
  msel_res$AIC <- est_paras_ID$AIC
  msel_res$Delta_AIC <- est_paras_ID$AIC - min(est_paras_ID$AIC)
  msel_res$AIC_corr <- est_paras_ID$AIC_corr
  msel_res$Delta_AIC_corr <- est_paras_ID$AIC_corr - min(est_paras_ID$AIC_corr)
  msel_res$p_value_UT <- est_paras_ID$p_value_UT
  msel_res$p_value_S <- est_paras_ID$p_value_S
  msel_res$p_value_joint <- est_paras_ID$p_value_joint
  write.csv(msel_res, file = paste0("experimental_data/model_selection/est_sum_", id, ".csv"))
}
