meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]

est_model <- c("no_SIM", "hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr", "het_zero_div", "het_div")
fit_m <- list("no_SIM"=1, "hom_wo_fitm"=1, "hom_fitm"=FALSE, "hom_fitm_unconstr"=c(FALSE,FALSE), "het_zero_div"=1, "het_div"=1)
rel_div_on <- list("no_SIM"=0, "hom_wo_fitm"=0, "hom_fitm"=0, "hom_fitm_unconstr"=0, "het_zero_div"=0, "het_div"=FALSE)
mod <- list("no_SIM"="no SIM", "hom_wo_fitm"="homogeneous", "hom_fitm"="homogeneous", "hom_fitm_unconstr"="homogeneous", "het_zero_div"="heterogeneous", "het_div"="heterogeneous")
est_p <- c("mu_UT","fitm_UT","mu_S","fitm_S", "fitm_ratio", "M", "mu_off", "S", "mu_on", "f_on", "rel_div_on", "mu_inc")
pq <- c()
for (p in est_p) {
  for (q in c("MLE", "lower_bound", "upper_bound")) {
    pq <- append(pq, paste0(p, "_", q))
  }  
}
est_paras <- data.frame(matrix(ncol = (length(pq)+6), nrow = 0))
colnames(est_paras) <- c("ID", "model", "status", pq, c("LL", "AIC", "BIC"))

for (i in 1:length(meta_data$ID)) {
  mc_data <- read.table(paste0("experimental_data/raw_counts/", meta_data$ID[i], ".txt"), header = FALSE, sep = ",", fill = TRUE)
  mc_S <- read_counts(mc_data[2,])
  Nf_S <- read_counts(mc_data[3,])
  eff_S <- as.numeric(mc_data[4,1])
  mc_data_b <- read.table(paste0("experimental_data/raw_counts/", meta_data$baseline_ID[i], ".txt"), header = FALSE, sep = ",", fill = TRUE)
  mc_UT <- read_counts(mc_data_b[2,])
  Nf_UT <- read_counts(mc_data_b[3,])
  eff_UT <- as.numeric(mc_data_b[4,1])
  print(meta_data$ID[i])
  for (m in est_model) {
    res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, eff = c(eff_UT, eff_S), fit_m = fit_m[[m]], rel_div_on = rel_div_on[[m]], mod = mod[[m]]) 
    if (is.data.frame(res[[2]]) && res[[2]]$LL != -Inf) {
      if (mod[m] == "no SIM"){
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:2,4:6])), c(t(res[[1]][1,4:6])), rep(1,9), c(t(res[[1]][1,4:6])), rep(0,3), rep(NA,12), c(t(res[[2]][1,3:5])))
      }
      if (mod[m] == "homogeneous") {
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:6,4:6])), rep(NA,18), c(t(res[[2]][1,3:5])))
      }
      if (m == "het_div") {
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,15), c(t(res[[1]][9,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][4:8,4:6])), c(t(res[[2]][1,3:5])))
      }
      if (m == "het_zero_div") {
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,18), c(t(res[[1]][1:2,4:6])), rep(NA,6), c(t(res[[1]][3,4:6])), rep(NA,3), c(t(res[[2]][1,3:5])))
      }
    } else {
      est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "failed", rep(NA, 36), -Inf, -Inf, -Inf)
    }
  }
}


est_paras[,4:dim(est_paras)[2]][!is.na(est_paras[,4:dim(est_paras)[2]])] <- as.numeric(est_paras[,4:dim(est_paras)[2]][!is.na(est_paras[,4:dim(est_paras)[2]])])
#est_paras <- arrange(est_paras, ID)
write.csv(est_paras, file = "experimental_data/est_paras.csv")
