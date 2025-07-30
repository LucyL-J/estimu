library("ggplot2")

# Frenoy et al. norfloxacin: the homogeneous model with unconstrained mutant fitness (HOM2) is selected 
df_msel_Nor <- read.csv("experimental_data/model_selection/est_sum_Frenoy_Nor.csv")[,-1]
print(df_msel_Nor)
# Looking at this model selection output:
# The model with the lowest AIC_corr is HOM2, with the second lowest is HET2 (Delta_AIC_corr = 12.5)

# Reading the input data
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_LB.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])       # Mutant counts
Nf_UT <- mean(read_counts(mc_data[3,])) # Average final population size
eff_UT <- as.numeric(mc_data[4,1])      # Plating efficiency
# Norfloxacin condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_Nor.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])       # Mutant counts
Nf_S <- mean(read_counts(mc_data[3,])) # Average final population size
eff_S <- as.numeric(mc_data[4,1])      # Plating efficiency

# Using the homogeneous model without differential mutant fitness, HOM0, as in the origional study by Frenoy et al.
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), mod = "hom")
paras_HOM0 <- res[[1]]
print(paras_HOM0)
msel_HOM0 <- res[[2]]
print(msel_HOM0)

# Estimation under model HOM2 (by setting fit_m = c(FALSE, FALSE) the mutant fitness is inferred in both UT and S conditions separately)
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = c(FALSE, FALSE), mod = "hom")
# Data frame with the estimated parameters
paras_HOM2 <- res[[1]]
print(paras_HOM2)
# Data frame with the results relevant for model selection
msel_HOM2 <- res[[2]]
print(msel_HOM2)

# Estimation under model HET2 (by setting f_on = FALSE and rel_div_on = FALSE, both are inferred in the S condition)
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), f_on = FALSE, rel_div_on = FALSE, mod = "het")
paras_HET2 <- res[[1]]
print(paras_HET2)
msel_HET2 <- res[[2]]
print(msel_HET2)

# Model HOM2 estimates a mutant fitness advantage in the UT condition, which we don't account for in model HET2 -> poor model fit of HET2 in the UT condition
# We can use the estimates from HOM2 in the estimation of HET2
fit_m_UT <- paras_HOM2$MLE[2]
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = c(fit_m_UT, fit_m_UT), f_on = FALSE, rel_div_on = FALSE, mod = "het")
paras_HET2_fitm <- res[[1]]
print(paras_HET2_fitm)
msel_HET2_fitm <- res[[2]]
print(msel_HET2_fitm)

# Comparing AIC_corrs again
msel_HOM2$AIC_corr[1]      # HOM2
msel_HET2$AIC_corr[1]      # HET2
msel_HET2_fitm$AIC_corr[1] # HET2 with differential mutant fitness from HOM2 (UT condition)

# Estimation under the heterogeneous-response model, setting the fraction of on-cells to some exemplary values
# HET1 because relative division rate of on-cells is inferred
# Fraction = 64%, as estimated by Bulssico et al. for 10 ng/mL of ciprofloxacin
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = c(fit_m_UT, fit_m_UT), f_on = 0.64, rel_div_on = FALSE, mod = "het")
paras_HET1_fitm_64 <- res[[1]]
print(paras_HET1_fitm_64)
msel_HET1_fitm_64 <- res[[2]]
print(msel_HET1_fitm_64)
# Fraction = 5%, as estimated by Jaramillo-Riveri et al. for 3 ng/mL of ciprofloxacin
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = c(fit_m_UT, fit_m_UT), f_on = 0.05, rel_div_on = FALSE, mod = "het")
paras_HET1_fitm_5 <- res[[1]]
print(paras_HET1_fitm_5)
msel_HET1_fitm_5 <- res[[2]]
print(msel_HET1_fitm_5)

# Comparing estimated increase in mutation rate associated with the induction of the stress response
paras_HET2_fitm$MLE[8]    # Heterogeneous with f_on and rel_div_on both inferred
paras_HET1_fitm_64$MLE[8] # Heterogeneous with f_on=64% and rel_div_on inferred
paras_HET1_fitm_5$MLE[8]  # Heterogeneous with f_on=5% and rel_div_on inferred

# And comparing AIC_corrs again
msel_HOM2$AIC_corr[1]         # HOM2
msel_HET2$AIC_corr[1]         # HET2
msel_HET2_fitm$AIC_corr[1]    # HET2 with differential mutant fitness from HOM2 (UT condition) and f_on inferred
msel_HET1_fitm_64$AIC_corr[1] # HET2 with differential mutant fitness from HOM2 (UT condition) and f_on=64%
msel_HET1_fitm_5$AIC_corr[1]  # HET2 with differential mutant fitness from HOM2 (UT condition) and f_on=5%

# If we constrain the differential mutant fitness to be the same under UT and S conditions, the best homogeneous model is HOM1
# Estimation under model HOM1 (by setting fit_m = FALSE the mutant fitness is jointly inferred in UT and S conditions)
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = FALSE, mod = "hom")
paras_HOM1 <- res[[1]]
print(paras_HOM1)
msel_HOM1 <- res[[2]]
print(msel_HOM1)

# UT condition: calculate the mutant count distributions under HOM2, HOM1, HET2 and HET2 with mutant fitness as estimated by HOM2
p_mc_HOM2_UT <- pMudi(max(mc_UT), Nf_UT, paras_HOM2$MLE[1], plateff=eff_UT, fit_m=paras_HOM2$MLE[2]) * length(mc_UT)
p_mc_HOM1_UT <- pMudi(max(mc_UT), Nf_UT, paras_HOM1$MLE[1], plateff=eff_UT, fit_m=paras_HOM1$MLE[2]) * length(mc_UT)
p_mc_HET2_UT <- pMudi(max(mc_UT), Nf_UT, paras_HET2$MLE[1], plateff=eff_UT) * length(mc_UT)
p_mc_HET2_fitm_UT <- pMudi(max(mc_UT), Nf_UT, paras_HET2_fitm$MLE[1], plateff=eff_UT, fit_m=paras_HET2_fitm$MLE[2]) * length(mc_UT)
# Estimated mutation rate and mutant fitness under HOM2, and estimated mutation rate and mutant fitness under HOM1
print(c(paras_HOM2$MLE[1], paras_HOM2$MLE[2], paras_HOM1$MLE[1], paras_HOM1$MLE[2]))
# Estimated mutation rate off-cells and fixed mutant fitness under HET2, and estimated mutation rate off-cells and fixed mutant fitness under HET2 with mutant fitness from HOM2
print(c(paras_HET2$MLE[1], paras_HET2$MLE[2], paras_HET2_fitm$MLE[1], paras_HET2_fitm$MLE[2]))
p_mc_UT <- ggplot() + geom_histogram(aes(mc_UT, y=after_stat(density) * length(mc_UT)), fill="#1C458A", bins = 50) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HOM2_UT), color="#F0E442", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HOM1_UT), color="#E69F00", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HET2_UT), color="#CC79A7", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HET2_fitm_UT), color="#7570B3", linewidth = 1., linetype = "dashed") + 
  ggtitle("UT condition") + xlab("Number of colonies") + ylab("Number of plates") + ylim(-0.1, 7.1)
p_mc_UT

# S condition: 0.05 mug/mL norfloxacin: calculate the mutant count distributions under the same models as above
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_Nor.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])
Nf_S <- mean(read_counts(mc_data[3,]))
eff_S <- as.numeric(mc_data[4,1])
p_mc_HOM2_S <- pMudi(max(mc_S), Nf_S, paras_HOM2$MLE[3], plateff=eff_S, fit_m=paras_HOM2$MLE[4]) * length(mc_S)
p_mc_HOM1_S <- pMudi(max(mc_S), Nf_S, paras_HOM1$MLE[3], plateff=eff_S, fit_m=paras_HOM1$MLE[4]) * length(mc_S)
p_mc_HET2_S <- pMudi(max(mc_S), Nf_S, paras_HET2$MLE[1], plateff=eff_S, S=paras_HET2$MLE[4], f_on=paras_HET2$MLE[6], rel_div_on=paras_HET2$MLE[5]) * length(mc_S)
p_mc_HET2_fitm_S <- pMudi(max(mc_S), Nf_S, paras_HET2$MLE[1], plateff=eff_S, S=paras_HET2$MLE[4], f_on=paras_HET2$MLE[6], rel_div_on=paras_HET2$MLE[5], fit_m=paras_HET2_fitm$MLE[2]) * length(mc_S)
# Estimated mutation rate and mutant fitness under HOM2, and estimated mutation rate and mutant fitness under HOM1
print(c(paras_HOM2$MLE[3], paras_HOM2$MLE[4], paras_HOM1$MLE[3], paras_HOM2$MLE[4]))
# Estimated mutation rate off-cells and fixed mutant fitness, and estimated mutation-supply ratio, fraction on-cells and relative division rate on-cells under HET2
print(c(paras_HET2$MLE[1], paras_HET2$MLE[3], paras_HET2$MLE[4], paras_HET2$MLE[6], paras_HET2$MLE[5]))
# Estimated mutation rate off-cells and fixed mutant fitness, and estimated mutation-supply ratio, fraction on-cells and relative division rate on-cells under HET2 with mutant fitness cost from HOM2
print(c(paras_HET2_fitm$MLE[1], paras_HET2_fitm$MLE[3], paras_HET2_fitm$MLE[4], paras_HET2_fitm$MLE[6], paras_HET2_fitm$MLE[5]))
p_mc_s <- ggplot() + geom_histogram(aes(mc_S, y=..density.. *length(mc_S)), fill="#8EC44F", bins = 20) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_HOM2_S), color="#F0E442", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_HOM1_S), color="#E69F00", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_HET2_S), color="#CC79A7", linewidth = 1., linetype = "dashed") +
  geom_line(aes(x=0:max(mc_S), y=p_mc_HET2_fitm_S), color="#7570B3", linewidth = 1.) +
  xlab("Number of colonies") + ylab("Number of plates") + ggtitle("S condition (50 ng/mL norfloxacin)") + ylim(-0.5,34)
p_mc_s


# Mo et al. ciprofloxacin: the homogeneous model with unconstrained mutant fitness (HOM2) is selected 
df_msel_Cip <- read.csv("experimental_data/model_selection/est_sum_Mo_Cip_MG1655.csv")[,-1]
print(df_msel_Cip)
# Looking at this model selection output:
# The model with the lowest AIC_corr is HET0, with the second lowest is HOM0 (Delta_AIC_corr = 0.11)

# Reading the input data
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Mo_LB_MG1655.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])       # Mutant counts
Nf_UT <- mean(read_counts(mc_data[3,])) # Average final population size
eff_UT <- as.numeric(mc_data[4,1])      # Plating efficiency
# Cirpofloxacin condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Mo_Cip_MG1655.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])       # Mutant counts
Nf_S <- mean(read_counts(mc_data[3,])) # Average final population size
eff_S <- as.numeric(mc_data[4,1])      # Plating efficiency

# Estimation under model HOM1 (by setting fit_m = FALSE)
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = FALSE, mod = "hom")
# Data frame with the estimated parameters
paras_HOM1 <- res[[1]]
print(paras_HOM1)
# Data frame with the results relevant for model selection
msel_HOM1 <- res[[2]]
print(msel_HOM1)

# Estimation under model HET0 (by setting f_on = FALSE and rel_div_on = 0.)
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), f_on = FALSE, rel_div_on = 0., mod = "het")
# Data frame with the estimated parameters
paras_HET0 <- res[[1]]
print(paras_HET0)
# Data frame with the results relevant for model selection
msel_HET0 <- res[[2]]
print(msel_HET0)

# Estimation under model HETF0, setting the fraction of on-cells to some exemplary values
# Fraction = 64%, as estimated by Bulssico et al. for 10 ng/mL of ciprofloxacin
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), f_on = 0.64, rel_div_on = 0., mod = "het")
# Data frame with the estimated parameters
paras_HETF0_64 <- res[[1]]
print(paras_HETF0_64)
# Data frame with the results relevant for model selection
msel_HETF0_64 <- res[[2]]
print(msel_HETF0_64)
# Fraction = 5%, as estimated by Jaramillo-Riveri et al. for 3 ng/mL of ciprofloxacin
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), f_on = 0.05, rel_div_on = 0., mod = "het")
# Data frame with the estimated parameters
paras_HETF0_5 <- res[[1]]
print(paras_HETF0_5)
# Data frame with the results relevant for model selection
msel_HETF0_5 <- res[[2]]
print(msel_HETF0_5)

# Comparing estimated increase in mutation rate associated with the induction of the stress response
paras_HETF0_64$MLE[8] # Heterogeneous with f_on=64% and rel_div_on inferred
paras_HETF0_5$MLE[8]  # Heterogeneous with f_on=5% and rel_div_on inferred
# And comparing with the estimated increase in population-wide mutation rate 
paras_HETF0_64$MLE[9] # Heterogeneous with f_on=64% and rel_div_on inferred
paras_HETF0_5$MLE[9]  # Heterogeneous with f_on=5% and rel_div_on inferred
paras_HOM1$MLE[6]     # Homogeneous (without differential mutant fitness)

