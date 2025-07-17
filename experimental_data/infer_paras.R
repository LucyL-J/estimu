meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]

est_p <- c("mu_UT","fitm_UT","mu_S","fitm_S", "fitm_ratio", "M", "mu_off", "S", "rel_div_on", "f_on", "mu_on", "mu_inc")
pq <- c()
for (p in est_p) {
  for (q in c("MLE", "lower_bound", "upper_bound")) {
    pq <- append(pq, paste0(p, "_", q))
  }  
}
lk <- c()
for (k in c("joint", "UT", "S", "test")) {
  if (k == "joint") {
    for (l in c("LL", "p_value")) {
      lk <- append(lk, paste0(l, "_", k))
    }
  } else {
    for (l in c("LL", "p_value", "cutoff_count", "tail_prob")) {
      lk <- append(lk, paste0(l, "_", k))
    }
  }
}
est_paras <- data.frame(matrix(ncol = (length(pq)+length(lk))+6, nrow = 0))
colnames(est_paras) <- c("ID", "model", "status", pq, c("AIC", "AIC_corr"), lk, "calc_time")

est_model <- c("N0", "N1", "N2", "HOM0", "HOM1", "HOM2", "HET0", "HET2", "HETF0", "HET1")
fit_m <- list("N0"=1, "N1"=FALSE, "N2"=c(FALSE,FALSE), "HOM0"=1, "HOM1"=FALSE, "HOM2"=c(FALSE,FALSE), "HET0"=1, "HET2"=1, "HETF0"=1, "HET1"=1)
rel_div_on <- list("N0"=0, "N1"=0, "N2"=0, "HOM0"=0, "HOM1"=0, "HOM2"=0, "HET0"=0, "HET2"=FALSE, "HETF0"=0, "HET1"=FALSE)
mod <- list("N0"="null", "N1"="null", "N2"="null", "HOM0"="hom", "HOM1"="hom", "HOM2"="hom", "HET0"="het", "HET2"="het", "HETF0"="het", "HET1"="het")

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
    if(is.element(m, c("HETF0", "HET1"))){
      if(!is.na(as.numeric(meta_data$SOS_induction[i])) && is.na(as.logical(meta_data$SOS_induction[i]))){
        f_on <- as.numeric(meta_data$SOS_induction[i])
      } else {
        st <- FALSE
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "NA", rep(NA,36), Inf, Inf, -Inf, 0, rep(c(-Inf,0,0,1),2), rep(NA,4), 0)
      }
    }
    if(st){
      res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = fit_m[[m]], rel_div_on = rel_div_on[[m]], f_on = f_on, mod = mod[[m]]) 
      if (is.data.frame(res[[2]]) && res[[2]]$LL[1] != -Inf) {
        if (mod[m] == "null"){
          if (m == "N0") {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:2,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][3:4,4:6])), rep(1,3), c(t(res[[1]][1,4:6])), rep(0,3), rep(NA,12), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:3,6:9])), rep(NA,4), res[[2]][1,10])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:2,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][3:4,4:6])), rep(1,3), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:3,6:9])), rep(NA,4), res[[2]][1,10])
          }
        }
        if (mod[m] == "hom") {
          if (m == "HOM1"){
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:6,4:6])), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:3,6:9])), rep(NA,4), res[[2]][1,10])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", c(t(res[[1]][1:6,4:6])), rep(NA,18), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:4,6:9])), res[[2]][1,10]) 
          }
        }
        if (mod[m] == "het"){
          if (m == "HET0") {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,18), c(t(res[[1]][1,4:6])), c(t(res[[1]][4:5,4:6])), rep(NA,9), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:3,6:9])), rep(NA,4), res[[2]][1,10])
          } else {
            est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "success", rep(NA,15), c(t(res[[1]][9,4:6])), c(t(res[[1]][1,4:6])), c(t(res[[1]][4:8,4:6])), c(t(res[[2]][1,3:4])), c(t(res[[2]][1,6:7])), c(t(res[[2]][2:3,6:9])), rep(NA,4), res[[2]][1,10])
          }
        }
      } else {
        est_paras[nrow(est_paras) + 1,] <- c(meta_data$ID[i], m, "failed", rep(NA,36), Inf, Inf, -Inf, 0, rep(c(-Inf,0,0,1),2), rep(NA,4), 0)
      } 
    }
  }
}

for (i in 4:ncol(est_paras)) { est_paras[,i] <- as.numeric(est_paras[,i]) }
write.csv(est_paras, file = "experimental_data/est_paras.csv")