meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]

est_model <- c("no_SIM_wo_fitm", "no_SIM_fitm", "no_SIM_fitm_unconstr", "hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr", "het_zero_div", "het_div", "het_zero_div_fon", "het_div_fon")
fit_m <- list("no_SIM_wo_fitm"=1, "no_SIM_fitm"=FALSE, "no_SIM_fitm_unconstr"=c(FALSE,FALSE), "hom_wo_fitm"=1, "hom_fitm"=FALSE, "hom_fitm_unconstr"=c(FALSE,FALSE), "het_zero_div"=1, "het_div"=1, "het_zero_div_fon"=1, "het_div_fon"=1)
rel_div_on <- list("no_SIM_wo_fitm"=0, "no_SIM_fitm"=0, "no_SIM_fitm_unconstr"=0, "hom_wo_fitm"=0, "hom_fitm"=0, "hom_fitm_unconstr"=0, "het_zero_div"=0, "het_div"=FALSE, "het_zero_div_fon"=0, "het_div_fon"=FALSE)
mod <- list("no_SIM_wo_fitm"="null", "no_SIM_fitm"="null", "no_SIM_fitm_unconstr"="null", "hom_wo_fitm"="homogeneous", "hom_fitm"="homogeneous", "hom_fitm_unconstr"="homogeneous", "het_zero_div"="heterogeneous", "het_div"="heterogeneous", "het_zero_div_fon"="heterogeneous", "het_div_fon"="heterogeneous")
est_p <- c("mu_UT","fitm_UT","mu_S","fitm_S", "fitm_ratio", "M", "mu_off", "S", "mu_on", "f_on", "rel_div_on", "mu_inc")
pq <- c()
for (p in est_p) {
  for (q in c("MLE", "lower_bound", "upper_bound")) {
    pq <- append(pq, paste0(p, "_", q))
  }  
}
lcond <- c()
for (cond in c("joint", "UT", "S", "test")) {
  for (l in c("LL", "p_value")) {
    lcond <- append(lcond, paste0(l, "_", cond))
  }
}
est_paras <- data.frame(matrix(ncol = (length(pq)+5+length(lcond))+1, nrow = 0))
colnames(est_paras) <- c("ID", "model", "status", pq, c("AIC_joint", "BIC_joint"), lcond, "time")

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
    st <- TRUE
    f_on <- FALSE
    if(is.element(m, c("het_zero_div_fon", "het_div_fon"))){
      if(!is.na(as.numeric(meta_data$SOS_induction[i])) && is.na(as.logical(meta_data$SOS_induction[i]))){
        f_on <- as.numeric(meta_data$SOS_induction[i])
      } else {
        st <- FALSE
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "NA", rep(NA,36), -Inf, -Inf, rep(c(-Inf,0),3), rep(NA,2), 0)
      }
    }
    if(st){
      res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = fit_m[[m]], rel_div_on = rel_div_on[[m]], f_on = f_on, mod = mod[[m]]) 
      if (is.data.frame(res[[2]]) && res[[2]]$LL[1] != -Inf) {
        if (mod[m] == "null"){
          if (m == "no_SIM_wo_fitm") {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:2,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][3:4,4:6])), rep(1,3), c(t(res[[1]][1,4:6])), rep(0,3), rep(NA,12), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:3,5:6])), rep(NA,2), res[[2]][1,7])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:2,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][3:4,4:6])), rep(1,3), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:3,5:6])), rep(NA,2), res[[2]][1,7])
          }
        }
        if (mod[m] == "homogeneous") {
          if (m == "hom_fitm"){
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:6,4:6])), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:3,5:6])), rep(NA,2), res[[2]][1,7])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:6,4:6])), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:4,5:6])), res[[2]][1,7]) 
          }
        }
        if (mod[m] == "heterogeneous"){
          if (m == "het_zero_div") {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,18), c(t(res[[1]][1,4:6])), c(t(res[[1]][4,4:6])), rep(NA,6), rep(0,3), rep(NA,3), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:3,5:6])), rep(NA,2), res[[2]][1,7])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,15), c(t(res[[1]][9,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][4:8,4:6])), c(t(res[[2]][1,3:4])), c(t(res[[2]][1:3,5:6])), rep(NA,2), res[[2]][1,7])
          }
        }
      } else {
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "failed", rep(NA,36), Inf, Inf, rep(c(-Inf,0),3), rep(NA,2), res[[2]][1,7])
      } 
    }
  }
}

est_paras[,4:dim(est_paras)[2]][!is.na(est_paras[,4:dim(est_paras)[2]])] <- as.numeric(est_paras[,4:dim(est_paras)[2]][!is.na(est_paras[,4:dim(est_paras)[2]])])
write.csv(est_paras, file = "experimental_data/est_paras.csv")
