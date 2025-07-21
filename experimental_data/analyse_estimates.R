library("plyr")
library("ggplot2")
library("viridisLite")
library("ggpubr")
library("lme4")
library("latex2exp")

# Read data frames
antibiotic_classes <- read.csv("experimental_data/antibiotic_classes.csv")[,-1]
meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]
sum_data <- read.csv("experimental_data/sum_data.csv")[,-1]
est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
est_sum <- read.csv("experimental_data/est_sum.csv")[,-1]

# Add the target (as we categorise it) to the meta data
meta_data$target <- mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$target_group)
# Use pooled experiments only (and not the individual replicates)
meta_data <- subset(meta_data, replicate == 0)
# The master data frame with all information
df <- merge(meta_data, sum_data, by = "ID")
df <- merge(df, est_sum, by = "ID")
# Make certain columns into factors, for plotting 
df$antibiotic <- factor(df$antibiotic, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
df$target <- factor(df$target, levels = unique(antibiotic_classes$target_group), ordered = TRUE)
antibiotic_classes$antibiotic_abbr <- factor(antibiotic_classes$antibiotic_abbr, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
antibiotic_classes$target_group <- factor(antibiotic_classes$target_group, levels = unique(antibiotic_classes$target_group), ordered = TRUE)
df$hom_by_AIC <- factor(df$hom_by_AIC, levels = c("N0", "N1", "N2", "HOM0", "HOM1", "HOM2"))
df$hom_by_AIC_corr <- factor(df$hom_by_AIC_corr, levels = c("N0", "N1", "N2", "HOM0", "HOM1", "HOM2"))

# Color-coding antibiotics used in the studies (sorted by target group)
antibiotic_classes$color <- turbo(length(antibiotic_classes$antibiotic_abbr))
# How many experiments use a certain antibiotic?
v <- numeric(length(antibiotic_classes$antibiotic_abbr))
for (i in 1:length(v)) { v[i] <- length(subset(meta_data, antibiotic == antibiotic_classes$antibiotic_abbr[i])$antibiotic) }
antibiotic_classes$prevalence <- v
p_antibiotic <- ggplot(data = antibiotic_classes, aes(x=target_group, y=prevalence, fill=antibiotic_abbr)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = antibiotic_classes$color, name = "Antimicrobial (abbr)") + 
  xlab("Grouped by target") + ylab("Number of experiments") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.8))
p_antibiotic

# Sort by bacterial species, for plotting
df <- arrange(df, species, target)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)

# Goodness-of-fit test for the UT condition (standard model without or with differential mutant fitness)
# Print the experiments that fail both GoF tests (standard without/with differential mutant fitness, p<0.05)
print(as.character(subset(df, p_value_UT_max < 0.05)$ID))
# Restrict to experiments that pass at least one GoF test (p>=0.05, all experiments pass as of 17/07/2025)
df <- subset(df, p_value_UT_max >= 0.05)


# Main figure Fig1: Estimated increase in population-wide mutation rate for all experiments (in various versions)

# Fig1 (i) Model used in the inference: homogeneous-response without differential mutant fitness (HOM0), 
# for purpose of comparison with previous results, as this model would have been used in previous studies
p_M_HOM0 <- ggplot(data = df, aes(x=ID, y=M_HOM0.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic), shape = 17) +
  geom_errorbar(aes(ymin=M_HOM0.2, ymax=M_HOM0.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + #guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_HOM0
# Print experiments for which an increase in mutation rate is estimated (model = homogeneous-response without differential mutant fitness)
print(length(subset(df, M_HOM0.1 > 1)$ID))


# Do the width of the confidence intervals of the estimated change in population-wide mutation rate M depend on the experimental desgin, 
# i.e. plated fraction and number of parallel cultures
# Print: how many experiments use a plated fraction smaller than 1
print(length(subset(df, plated_fraction < 1)$ID))
df$width_CI <- (df$M_HOM0.3-df$M_HOM0.2)/df$M_HOM0.1
# Fitting a linear model with the log-width of the confidence interval of M as predicted variable, and the log-plated fraction and log-number of parallel cultures as predictor
lm <- lm(log10(width_CI) ~ log10(plated_fraction) + log10(n_cultures), data = df)
summary(lm)
p_CI <- ggplot(data = df, aes(x=(plated_fraction), y=(width_CI))) + #+log10(n_cultures)
  geom_point(aes(color=log10(n_cultures))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans="log10") +
  labs(x="Plated fraction", y="Normalised width of 95% CI around MLE estimate", color=TeX("$log_{10}(c_s)$"))
p_CI
# Correlation between the plated fraction and the number of parallel cultures
cor.test(df$plated_fraction, df$n_cultures, method = "kendall")


# Fig1 (ii) Model used in the inference: the homogeneous-response/null model with the lowest AIC (N0, N1, N2, HOM0, HOM1 or HOM2)
p_M_AIC <- ggplot(data = df, aes(x=ID, y=M_AIC.1, group=antibiotic, shape = hom_by_AIC)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +  scale_shape_manual(values = c(2,5,0,17,18,15), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_AIC
# Print experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC)
print(length(subset(df, M_AIC.1 > 1)$ID))

# Fig1 (iii) Model used in the inference: the homogeneous-response/null model with the lowest AIC_corr
# AIC corrected for the sample size n: AIC_corr := 2p n/(n-p-1) - 2ln(L)
p_M_AIC_corr <- ggplot(data = df, aes(x=ID, y=M_AIC_corr.1, group=antibiotic, shape = hom_by_AIC_corr)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_AIC_corr
# Print experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC_corr)
print(length(subset(df, M_AIC_corr.1 > 1)$ID))

# Null hypothesis H0: the data observed in the S condition could have been generated under the parameters inferred from the UT condition
# We test for this using a goodness-of-fit test

# Print: experiments for which we reject the null hypothesis H0 (criterion: p<0.05 in GoF)
print(length(subset(df, p_value_test_min < 0.05)$ID))

# Comparison: experiments with an estimated change in mutation rate (M!=1), but we do not reject H0
setdiff(subset(df, M_AIC.1 != 1)$ID, subset(df, p_value_test_min >= 0.05)$ID)      # AIC
setdiff(subset(df, M_AIC_corr.1 != 1)$ID, subset(df, p_value_test_min >= 0.05)$ID) # AIC_corr
# and the other way round
setdiff(subset(df, p_value_test_min >= 0.05)$ID, subset(df, M_AIC.1 != 1)$ID)      # AIC
setdiff(subset(df, p_value_test_min >= 0.05)$ID, subset(df, M_AIC_corr.1 != 1)$ID) # AIC_corr

# Does it depend on the experimental design, whether we accept/reject H0?
df_GoF <- df #subset(df, p_value_test_min > 0)
# Fitting a linear model with the p-value as predicted variable, and log-plated fraction and log-number of cultures as predictors
lm_p <- lm(p_value_test_min ~ log10(plated_fraction) + log10(n_cultures), data = df_GoF)
summary(lm_p)
p_GoF <- ggplot(data = df_GoF, aes(x=(n_cultures), y=(p_value_test_min))) + #+log10(n_cultures)
  geom_point(aes(color=log10(plated_fraction))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + #scale_y_continuous(trans = "log10") +
  labs(x="Number of parallel cultures", y="p value GoF test (S generated under UT parameters)", color=TeX("$log_{10}(E)$"))
p_GoF

# Restriction to experiments for which we reject H0 (setting the estimated change in mutation rate to M=1 if H0 is accepted)
df_strict <- df
print(as.character(subset(df_strict, p_value_test_min >= 0.05)$ID))
df_strict$M_AIC.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.3[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.3[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_HOM0.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_HOM0.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_HOM0.3[df_strict$p_value_test_min >= 0.05] <- 1.
levels(df_strict$hom_by_AIC) <- c(levels(df_strict$hom_by_AIC), "Failed GoF test")
df_strict$hom_by_AIC[df_strict$p_value_test_min >= 0.05] <- "Failed GoF test"
levels(df_strict$hom_by_AIC_corr) <- c(levels(df_strict$hom_by_AIC_corr), "Failed GoF test")
df_strict$hom_by_AIC_corr[df_strict$p_value_test_min >= 0.05] <- "Failed GoF test"
df_strict$H0 <- "HOM0"
df_strict$H0[df_strict$p_value_test_min >= 0.05] <- "Failed GoF test"

# Evaluate the model fits using goodness-of-fit tests
# Print: experiments for which the model with the lowest AIC is a poor fit to the data (none as of 21/07/2025)
print(as.character(subset(df_strict, p_value_UT_hom_AIC < 0.05)$ID)) # AIC in the UT condition
print(as.character(subset(df_strict, p_value_S_hom_AIC < 0.05)$ID))  # AIC in the S condition
print(as.character(subset(df_strict, p_value_UT_hom_AIC_corr < 0.05)$ID)) # AIC_corr in the UT condition
print(as.character(subset(df_strict, p_value_S_hom_AIC_corr < 0.05)$ID))  # AIC_corr in the S condition

# Are p-values of the GoF tests different for the UT and S conditions?
# My intuition here is that in the S condition it should be more likely that something completely different is going on, 
# which cannot be explained well by any of the models, not even the model with the lowest AIC_corr
df_p_value_UT_S <- data.frame(ID = rep(df$ID, 2))
df_p_value_UT_S$condition <- rep(c("UT", "S"), each = length(df$ID))
df_p_value_UT_S$p_values <- c(df$p_value_UT_hom_AIC_corr, df$p_value_S_hom_AIC_corr) # Take het into account!!!
wilcox.test(p_values ~ condition, data = df_p_value_UT_S)
median(subset(df_p_value_UT_S, condition == "UT")$p_values)
median(subset(df_p_value_UT_S, condition == "S")$p_values)
p_M_DNA <- ggplot(data = df_p_value_UT_S, aes(x=condition, y=p_values)) + 
  geom_boxplot(aes(fill=condition), show.legend = FALSE, outlier.shape = NA) + 
  ylab("p-value goodness of fit test") + xlab("Treatment condition") + ggtitle("Hom model with lowest AIC_corr") + 
  stat_compare_means(label.y = 1) + theme(legend.position = "right")
p_M_DNA

# Same as Fig1 version (i), but marking experiments with H0 not rejected, or poor model fits
p_M_HOM0_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_HOM0.1, group=antibiotic, shape = factor(H0))) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_HOM0.2, ymax=M_HOM0.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(8,17), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_HOM0_GoF

# Same as Fig1 version (ii) but marking experiments with H0 not rejected, or poor model fits
p_M_AIC_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_AIC.1, group=antibiotic, shape = hom_by_AIC)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15,8), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_AIC_GoF

# Same as Fig1 version (iii) but marking experiments with H0 not rejected, or poor model fits
p_M_AIC_corr_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_AIC_corr.1, group=antibiotic, shape = hom_by_AIC_corr)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15,8), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_AIC_corr_GoF


# Define 'stress-induced mutagenesis is detected' (from least to most strict)
# 1. The model with the lowest AIC estimates an increase in population-wide mutation rate M>1 (the MLE, not the bounds of the CI)
df$SIM_AIC <- df$M_AIC.1 > 1
print(length(subset(df, SIM_AIC == TRUE)$ID))
# 2. The model with the lowest AIC_corr estimates an increase M>1
df$SIM_AIC_corr <- df$M_AIC_corr.1 > 1
print(length(subset(df, SIM_AIC_corr == TRUE)$ID))
# 3. The model with the lowest AIC estimates an increase M>1, and the null hypothesis H0, that S data was generated under UT parameters, is rejected
df_strict$SIM_AIC <- df_strict$M_AIC.1 > 1
print(length(subset(df_strict, SIM_AIC == TRUE)$ID))
# 4. The model with the lowest AIC_corr estimates an increase M>1, and H0 is rejected
df_strict$SIM_AIC_corr <- df_strict$M_AIC_corr.1 > 1
print(length(subset(df_strict, SIM_AIC_corr == TRUE)$ID))

# Comparing to for how many experiments SIM would have been detected solely on the basis of estimating an increase M>1
print(length(subset(df, M_HOM0.1 > 1)$ID))

# We restrict the analysis about the target to experiments using E. coli MG1655 and TD2158 (no mutant strains)
df_Ecoli <- subset(subset(df, species == "E. coli"), is.element(strain, c("MG1655", "TD2158")))
df_Ecoli_strict <- subset(subset(df_strict, species == "E. coli"), is.element(strain, c("MG1655", "TD2158")))

# How many E. coli experiments use DNA/DNA gyrase-targeting antibiotics?
length(subset(df_Ecoli, is.element(target, c("DNA", "DNA gyrase")))$ID)
# And how many use ribosome-targeting ones?
length(subset(df_Ecoli, target=="Ribosome")$ID)

# Testing for normality of the estimated increase M
shapiro.test(df_Ecoli$M_HOM0.1)                                                     # All E. coli experiments
shapiro.test(subset(df_Ecoli, is.element(target, c("DNA", "DNA gyrase")))$M_HOM0.1) # E. coli experiments using a DNA/DNA gyrase-targeting antibiotic
shapiro.test(subset(df_Ecoli, target == "Ribosome")$M_HOM0.1)                       # E. coli experiments using a ribosome-targeting antibiotic

# Kruskal-Wallis test -> Does the estimated M depend on the antibiotic target?
kruskal.test(M_HOM0.1 ~ target, data = df_Ecoli)
# More specifically, is M significantly different for DNA/DNA gyrase vs. ribosome-targeting antibiotics?
df_KW <- subset(df_Ecoli, is.element(target, c("DNA", "DNA gyrase", "Ribosome")))
df_KW$group <- character(length(df_KW$target))
df_KW$group[is.element(df_KW$target, c("DNA", "DNA gyrase"))] <- "DNA/DNA gyrase"
df_KW$group[df_KW$target == "Ribosome"] <- "Ribosome"
wilcox.test(M_HOM0.1 ~ group, data = df_KW)
median(subset(df_KW, group == "DNA/DNA gyrase")$M_HOM0.1)
median(subset(df_KW, group == "Ribosome")$M_HOM0.1)
# We also plot experiments for which SIM is detected here (defined as model with lowest AIC_corr estimates M>1 and H0 is rejected)
df_KW$SIM <- subset(df_Ecoli_strict, is.element(target, c("DNA", "DNA gyrase", "Ribosome")))$SIM_AIC_corr
# For how many experiments do we detect SIM, depending on the target?
sum(subset(df_KW, group == "DNA/DNA gyrase")$SIM) # DNA/DNA gyrase-targeting 
sum(subset(df_KW, group == "Ribosome")$SIM)       # Ribosome-targeting
p_M_DNA <- ggplot(data = df_KW, aes(x=group, y=M_HOM0.1)) + geom_boxplot(aes(fill=group), show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(aes(color = SIM), width = 0.25, alpha = 0.8) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  scale_fill_manual(values = c("DNA/DNA gyrase" = "#4E6ADB", "Ribosome" = "#FE9B2D")) + 
  theme(plot.margin = margin(0.5,0.5,2.5,0.5, "cm")) + 
  ylab("Estimated fold change in population-wide mutation rate") + xlab("Antimicrobial target") +
  stat_compare_means(label.y = 300) + theme(legend.position = "right")
p_M_DNA
print(chisq.test(df_KW$group, df_KW$SIM_AIC_corr))

# What else does detection of SIM depend on?
# For our generalised linear mixed model, we again define SIM as 'model with the lowest AIC_corr estimates M>1 and H0 is rejected'
df_glmm <- df_strict
# The antibiotic dose in units of the MIC is a predictor variable, so we consider only experiments for which the MIC was measured
df_glmm <- subset(df_glmm, !is.na(of_MIC))
length(df_strict$ID) # All experiments
length(df_glmm$ID)   # Experiments for which the MIC was measured
# Grouping the target into 'DNA/DNA gyrase', 'Ribosome' and 'Other'
df_glmm$group <- character(length(df_glmm$target))
df_glmm$group[is.element(df_glmm$target, c("DNA", "DNA gyrase"))] <- "DNA/DNA gyrase"
df_glmm$group[df_glmm$target == "Ribosome"] <- "Ribosome"
df_glmm$group[!is.element(df_glmm$target, c("DNA", "DNA gyrase", "Ribosome"))] <- "Other"
df_glmm$group <- relevel(factor(df_glmm$group), ref = "Ribosome")
# Fitting a linear model with detection of SIM as a predicted variable
# and the antibiotic dose, the target group and the log-plated fraction as fixed effects
# and the baseline ID as random effect
glmm <- glmer(SIM_AIC_corr ~ of_MIC + group + log10(plated_fraction) + (1|baseline_ID), data = df_glmm, family=binomial)
summary(glmm)

# Fitting a linear mixed model with the p-value of testing H0 as a predicted variable 
# and the antibiotic dose, the target group and the log-plated fraction as fixed effects
# and the baseline ID as random effect
lmer_p <- lmer(p_value_test_min ~ of_MIC + group + log10(plated_fraction) + (1|baseline_ID), data = df_glmm)
summary(lmer_p)

# Experiments with E. coli for which SIM was detected (model with the lowest AIC_corr estimates M>1 and H0 is rejected)
df_SIM <- subset(df_Ecoli_strict, SIM_AIC_corr == TRUE)
length(df_SIM$ID)                       # Number of E. coli experiments with SIM
length(subset(df_Ecoli, M_HOM0.1>1)$ID) # Number of E. coli experiments with M>1 estimated by model HOM0 (previous studies)

# Estimating the mutation-supply ratio S (heterogeneous-response models) for experiments with detected SIM
# No model selection between homogeneous and heterogeneous models yet!

# Evaluate the heterogeneous model fits using goodness-of-fit tests
# Print: E. coli experiments for which the heterogeneous model with the lowest AIC is a poor fit to the data
print(as.character(subset(df_SIM, p_value_UT_het_AIC < 0.05)$ID)) # AIC in the UT condition
print(as.character(subset(df_SIM, p_value_S_het_AIC < 0.05)$ID))  # AIC in the S condition
print(as.character(subset(df_SIM, p_value_UT_het_AIC_corr < 0.05)$ID)) # AIC_corr in the UT condition
print(as.character(subset(df_SIM, p_value_S_het_AIC_corr < 0.05)$ID))  # AIC_corr in the S condition
# Print: all experiments for which the heterogeneous model with the lowest AIC is a poor fit to the data
print(as.character(subset(df, p_value_UT_het_AIC < 0.05)$ID)) # AIC in the UT condition
print(as.character(subset(df, p_value_S_het_AIC < 0.05)$ID))  # AIC in the S condition
print(as.character(subset(df, p_value_UT_het_AIC_corr < 0.05)$ID)) # AIC_corr in the UT condition
print(as.character(subset(df, p_value_S_het_AIC_corr < 0.05)$ID))  # AIC_corr in the S condition

# Generate separate data frames to mark experiments for which the model fits under UT/S conditions are poor
df_SIM_GoF <- subset(df_SIM, p_value_UT_het_AIC < 0.05 | p_value_S_het_AIC < 0.05)
levels(df_SIM_GoF$het_by_AIC) <- c(levels(df_SIM_GoF$het_by_AIC), "Poor fit UT cond", "Poor fit S cond")
df_SIM_GoF$het_by_AIC[df_SIM_GoF$p_value_UT_het_AIC < 0.05] <- "Poor fit UT cond"
df_SIM_GoF$het_by_AIC[df_SIM_GoF$p_value_S_het_AIC < 0.05] <- "Poor fit S cond"
df_SIM_GoF_corr <- subset(df_SIM, p_value_UT_het_AIC_corr < 0.05 | p_value_S_het_AIC_corr < 0.05)
levels(df_SIM_GoF_corr$het_by_AIC_corr) <- c(levels(df_SIM_GoF_corr$het_by_AIC_corr), "Poor fit UT cond", "Poor fit S cond")
df_SIM_GoF_corr$het_by_AIC_corr[df_SIM_GoF_corr$p_value_UT_het_AIC_corr < 0.05] <- "Poor fit UT cond"
df_SIM_GoF_corr$het_by_AIC_corr[df_SIM_GoF_corr$p_value_S_het_AIC_corr < 0.05] <- "Poor fit S cond"

p_S_SIM_all <- ggplot(data = df_SIM, aes(x=ID, y=S_AIC_corr.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic, shape = het_by_AIC_corr), size=2) + scale_shape_manual(values = c(1,16,2,4,3), name = "Selected model") +
  geom_errorbar(aes(ymin=S_AIC_corr.2, ymax=S_AIC_corr.3, color=antibiotic)) + geom_hline(yintercept = 1) +
  scale_colour_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Mutation-supply ratio") + xlab("Experiment ID") + geom_point(data = df_SIM_GoF_corr, aes(x=ID, y=S_AIC_corr.1, group=antibiotic, shape=het_by_AIC_corr), size = 3)
p_S_SIM_all

# Including model selection between selected homogeneous/heterogeneous-response models (criterion = lowest AIC_corr)
print(df_SIM$by_AIC_corr)
df_SIM <- arrange(df_SIM, target)
df_SIM$ID <- factor(df_SIM$ID, levels = unique(df_SIM$ID), ordered = TRUE)

# Estimated increase M for E. coli experiments with SIM
# Estimates for which a heterogeneous is selected over a homogeneous response overall are shown in dashed lines
df_SIM$homhet <- character(length(df_SIM$by_AIC_corr))
df_SIM$homhet[is.element(df_SIM$by_AIC_corr, c("HOM0", "HOM1", "HOM2"))] <- "Homogeneous"
df_SIM$homhet[is.element(df_SIM$by_AIC_corr, c("N0", "HET0", "HET2"))] <- "Heterogeneous"
df_SIM$homhet_sel <- df_SIM$hom_by_AIC_corr == df_SIM$by_AIC_corr
p_M_SIM <- ggplot(data = df_SIM, aes(x=ID, y=M_AIC_corr.1, group=target)) + 
  geom_point(aes(color=target, shape = hom_by_AIC_corr), size=2) + scale_shape_manual(values = c(17,18,15), name = "Selected model") +
  geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=target, linetype = homhet_sel)) + geom_hline(yintercept = 1) +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Overall selected") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate") + xlab("Experiment ID")
p_M_SIM

# Estimated increase S for E. coli experiments with SIM
# Estimates for which a homogeneous is selected over a heterogeneous response overall are shown in dashed lines
p_S_SIM <- ggplot(data = df_SIM, aes(x=ID, y=S_AIC_corr.1, group=target)) + 
  geom_point(aes(color=target, shape = het_by_AIC_corr), size=2) + scale_shape_manual(values = c(1,16,2), name = "Selected model") +
  geom_errorbar(aes(ymin=S_AIC_corr.2, ymax=S_AIC_corr.3, color=target, linetype = !homhet_sel)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Overall selected") +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") + #values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Mutation-supply ratio") + xlab("Experiment ID") #+ geom_point(data = df_SIM_GoF_corr, aes(x=ID, y=S_AIC_corr.1, group=antibiotic, shape=het_by_AIC_corr))
p_S_SIM

# Combining M und S plot with two y axes?
p_M_S_SIM <- ggplot(data = df_SIM, aes(x = ID)) +
  geom_point(aes(x = as.numeric(ID)-0.2, y = M_AIC_corr.1, color = target, shape = hom_by_AIC_corr), size = 2) +
  geom_errorbar(aes(x = as.numeric(ID)-0.2, ymin = M_AIC_corr.2, ymax = M_AIC_corr.3, color = target, linetype = homhet_sel), width = 0.5) +
  geom_point(aes(x = as.numeric(ID)+0.2, y = S_AIC_corr.1, color = target, shape = het_by_AIC_corr)) + 
  geom_errorbar(aes(x = as.numeric(ID)+0.2, ymin = S_AIC_corr.2, ymax = S_AIC_corr.3, color = target, linetype = !homhet_sel), width = 0.5) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Overall selected") +
  scale_shape_manual(values = c(17,18,15,1,16,2), name = "Selected model") +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  scale_y_continuous(trans="log10", "Increase population-wide mutation rate", sec.axis = sec_axis(~ ., name = "Mutation-supply ratio")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5, 0.5, 0.5, 0.5, "cm")) +
  xlab("Experiment ID") + scale_x_continuous(breaks = as.numeric(df_SIM$ID), labels = df_SIM$ID) #+ geom_hline(yintercept = 1)
p_M_S_SIM

# Model selection (criterion = lowest AIC_corr) between homogeneous- and heterogeneous-response models for E. coli experiments with SIM 
df_sel <- data.frame(antibiotic=rep(unique(df_SIM$antibiotic), each=2))
df_sel$target <- mapvalues(df_sel$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target_group))
df_sel$homhet <- factor(rep(c("Homogeneous","Heterogeneous"), length(unique(df_SIM$antibiotic))))
v <- numeric(length(df_sel$antibiotic))
for (i in 1:length(unique(df_SIM$antibiotic))) {
  v[2*i-1] <- sum(subset(df_SIM, antibiotic == df_sel$antibiotic[2*i])$homhet == "Homogeneous") 
  v[2*i] <- sum(subset(df_SIM, antibiotic == df_sel$antibiotic[2*i])$homhet == "Heterogeneous")
}
df_sel$prevalence <- v
p_sel_t <- ggplot(data = df_sel, aes(x=homhet, y=prevalence, fill=target)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  xlab("Overall selected") + ggtitle(criterion) + ylab("Number of experiments") + ggtitle("")#"E. coli experiments with detected SIM")
p_sel_t

# Difference in AIC_corr between homogeneous and heterogeneous-response models
p_Delta <- ggplot(data = df_SIM, aes(x=ID, y=Delta_AIC_corr, group=antibiotic)) + geom_point(aes(color=antibiotic, shape=by_AIC_corr)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  geom_hline(yintercept = 2, linetype = "dashed") + geom_hline(yintercept = -2, linetype = "dashed") + 
  ylab("Difference in AIC") + xlab("Experiment ID") + scale_shape_manual(values = c(1, 17,18,15), name = "Overall selected") 
p_Delta

# Case studies
i <- 20
id <- df_SIM$ID[i]
print(c(as.character(id), df_SIM$by_AIC[i]))
df_id <- read.csv(paste0("experimental_data/model_selection/est_sum_", id, ".csv"))[,-1]
print(as.character(subset(df_SIM, by_AIC == "HOM2")$ID))
print(as.character(subset(df_SIM, by_AIC == "HOM1")$ID))

# Experiments, for which a heterogeneous/homogeneous stress response is selected
df_het <- subset(df_SIM, is.element(by_AIC, c("het")))
df_hom <- subset(df_SIM, is.element(by_AIC, c("hom")))

# Frenoy et al. 2018 Norfloxacin
df_Nor <- subset(est_paras, ID == "Frenoy_Nor")
print(subset(est_sum, ID == "Frenoy_Nor"))
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_LB.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])
Nf_UT <- mean(read_counts(mc_data[3,]))
eff_UT <- as.numeric(mc_data[4,1])
p_mc_hom_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[6], plateff=eff_UT, fit_m=df_Nor$fitm_UT_MLE[6]) * length(mc_UT)
p_mc_hom_constr_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[4], plateff=eff_UT) * length(mc_UT)
p_mc_het_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_off_MLE[8], plateff=eff_UT) * length(mc_UT)
print(c(df_Nor$mu_UT_MLE[6], df_Nor$fitm_UT_MLE[6], df_Nor$mu_UT_MLE[4]))
print(df_Nor$mu_off_MLE[8])
p_mc_UT <- ggplot() + geom_histogram(aes(mc_UT, y=after_stat(density) * length(mc_UT)), fill="#1C458A", bins = 50) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_UT), color="#8F3F8C", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_constr_UT), color="#E0B0FF", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_het_UT), color="#FFD300", linewidth = 1., linetype = "dashed") + 
  ggtitle("Untreated experiments") + xlab("Number of colonies") + ylab("Number of plates") + ylim(-0.1, 7.1)
p_mc_UT
# Stressful condition: 0.05 mug/mL Norfloxacin
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_Nor.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])
Nf_S <- mean(read_counts(mc_data[3,]))
eff_S <- as.numeric(mc_data[4,1])
p_mc_hom_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[6], plateff=eff_S, fit_m=df_Nor$fitm_S_MLE[6]) * length(mc_S)
p_mc_hom_constr_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[4], plateff=eff_S) * length(mc_S)
p_mc_het_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_off_MLE[8], plateff=eff_S, S=df_Nor$S_MLE[8], f_on=df_Nor$f_on_MLE[8], rel_div_on=df_Nor$rel_div_on_MLE[8]) * length(mc_S)
print(c(df_Nor$mu_S_MLE[6], df_Nor$fitm_S_MLE[6], df_Nor$mu_S_MLE[4]))
print(c(df_Nor$mu_off_MLE[8], df_Nor$S_MLE[8], df_Nor$f_on_MLE[8], df_Nor$rel_div_on_MLE[8]))
p_mc_s <- ggplot() + geom_histogram(aes(mc_S, y=..density.. *length(mc_S)), fill="#8EC44F", bins = 20) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_S), color="#8F3F8C", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_constr_S), color="#E0B0FF", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_het_S), color="#FFD300", linewidth = 1., linetype = "dashed") +
  xlab("Number of colonies") + ylab("Number of plates") + ggtitle("Experiments with norfloxacin") + ylim(-0.5,34)
p_mc_s

# Dose-dependence
df_Pribis <- read.delim("experimental_data/Pribis_Fig2D.txt", header = TRUE, sep = '\t', comment.char="#")
df_c <- subset(df_Ecoli, antibiotic == "Cip")
p_c <- ggplot() + 
  geom_point(data = df_c, aes(x=concentration, y=M_HOM0.1)) +
  geom_errorbar(data = df_c, aes(x=concentration, ymin=M_HOM0.2, ymax=M_HOM0.3)) +
  geom_hline(yintercept = 1) + scale_y_continuous(trans="log10") +
  geom_point(data = df_Pribis, aes(x=cipro_concentration/1000, y=fold_induction_mutation_rate), color='blue') +
  geom_errorbar(data = df_Pribis, aes(x=cipro_concentration/1000, ymin=fold_min, ymax=fold_max), color='blue') +
  ylab("Fold-change population-wide mutation rate") + ggtitle("Cip")
p_c
