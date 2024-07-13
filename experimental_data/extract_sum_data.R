meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]

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

# Number of parallel cultures (_b respective baseline experiment)
n <- numeric(length(meta_data$ID))
n_b <- numeric(length(meta_data$ID))
# Plated fraction 
eff <- numeric(length(meta_data$ID))
eff_b <- numeric(length(meta_data$ID))
# Mean number and coefficient of variation of the final population size
mean_N <- numeric(length(meta_data$ID))
cv_N <- numeric(length(meta_data$ID))
mean_N_b <- numeric(length(meta_data$ID))
cv_N_b <- numeric(length(meta_data$ID))

for (i in 1:length(meta_data$ID)) {
  mc_data <- read.table(paste0("experimental_data/raw_counts/", meta_data$ID[i], ".txt"), header = FALSE, sep = ",", fill = TRUE)
  mc_data_b <- read.table(paste0("experimental_data/raw_counts/", meta_data$baseline_ID[i], ".txt"), header = FALSE, sep = ",", fill = TRUE)
  mc <- read_counts(mc_data[2,])
  Nf <- read_counts(mc_data[3,])
  mc_b <- read_counts(mc_data_b[2,])
  Nf_b <- read_counts(mc_data_b[3,])
  n[i] <- length(mc)
  n_b[i] <- length(mc_b)
  mean_N[i] <- mean(Nf)
  cv_N[i] <- (var(Nf))^0.5/mean_N[i]
  mean_N_b[i] <- mean(Nf_b)
  cv_N_b[i] <- (var(Nf_b))^0.5/mean_N_b[i]
  eff[i] <- as.numeric(mc_data[4,1])
  eff_b[i] <- as.numeric(mc_data_b[4,1])
}
n_tot <- n + n_b

sum_data <- data.frame(ID=meta_data$ID, n_cultures=n, n_cultures_baseline=n_b, n_cultures_tot=n_tot, plated_fraction=eff, plated_fraction_baseline=eff_b, mean_Nf=mean_N, CV_Nf=cv_N, mean_Nf_baseline=mean_N_b, CV_Nf_baseline=cv_N_b, replicate=meta_data$replicate, pooled = meta_data$pooled)
write.csv(sum_data, file = "experimental_data/sum_data.csv")