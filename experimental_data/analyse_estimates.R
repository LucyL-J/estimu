library("FSA")
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

# Add the target and target group (as we categorise it) to the meta data
meta_data$target <- plyr::mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$target)
meta_data$target_group <- plyr::mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$target_group)
# Use pooled experiments only (and not the individual replicates)
meta_data <- subset(meta_data, replicate == 0)
# The master data frame with all information
df <- merge(meta_data, sum_data, by = "ID")
df <- merge(df, est_sum, by = "ID")
# Make certain columns into factors, for plotting 
antibiotic_classes <- arrange(antibiotic_classes, target_group, target, antibiotic)
df$antibiotic <- factor(df$antibiotic, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
df$target <- factor(df$target, levels = unique(antibiotic_classes$target), ordered = TRUE)
antibiotic_classes$antibiotic_abbr <- factor(antibiotic_classes$antibiotic_abbr, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
antibiotic_classes$target <- factor(antibiotic_classes$target, levels = unique(antibiotic_classes$target), ordered = TRUE)
#antibiotic_classes$target_group <- factor(antibiotic_classes$target_group, levels = unique(antibiotic_classes$target_group), ordered = TRUE)
df$null_hom_by_AIC <- factor(df$null_hom_by_AIC, levels = c("N0", "N1", "N2", "HOM0", "HOM1", "HOM2"))
df$null_hom_by_AIC_corr <- factor(df$null_hom_by_AIC_corr, levels = c("N0", "N1", "N2", "HOM0", "HOM1", "HOM2"))
df$null_het_by_AIC <- factor(df$null_het_by_AIC, levels = c("N0", "HET0", "HET2"))
df$null_het_by_AIC_corr <- factor(df$null_het_by_AIC_corr, levels = c("N0", "HET0", "HET2"))

# Color-coding antibiotics used in the studies (sorted by target with a gap between groups)
g <- 3
all_colours <- turbo(length(antibiotic_classes$antibiotic)+1+2*g)
n_DNA <- length(subset(antibiotic_classes, target_group == "DNA/DNA gyrase")$antibiotic)
n_Other <- length(subset(antibiotic_classes, target_group == "Other")$antibiotic)
n_Ribosome <- length(subset(antibiotic_classes, target_group == "Ribosome")$antibiotic)
colours_w_gap <- c(all_colours[1:n_DNA], all_colours[(n_DNA+2+g):(n_DNA+n_Other+1+g)], all_colours[(n_DNA+n_Other+2+2*g):(n_DNA+n_Other+n_Ribosome+1+2*g)])
antibiotic_classes$color <- colours_w_gap
# How many experiments use a certain antibiotic?
v <- numeric(length(antibiotic_classes$antibiotic))
for (i in 1:length(v)) { v[i] <- length(subset(meta_data, antibiotic == antibiotic_classes$antibiotic_abbr[i])$antibiotic) }
antibiotic_classes$prevalence <- v
p_antibiotic <- ggplot(data = antibiotic_classes, aes(x=target, y=prevalence, fill=antibiotic_abbr)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = antibiotic_classes$color, name = "Antimicrobial (abbr)") + 
  xlab("Grouped by target") + ylab("Number of experiments") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.8))
p_antibiotic

# Sort by bacterial species, for plotting
df <- arrange(df, species, target_group, target, antibiotic)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)

# Goodness-of-fit test for the UT condition (standard model without or with differential mutant fitness)
# Print the experiments that fail both GoF tests (standard without/with differential mutant fitness, p<0.05)
print(as.character(subset(df, p_value_UT_max < 0.05)$ID))
# Restrict to experiments that pass at least one GoF test (p>=0.05, all experiments pass as of 24/07/2025)
df <- subset(df, p_value_UT_max >= 0.05)


# Results Figure 1: Estimated increase in population-wide mutation rate for all experiments (in various versions)

# Fig1 version (i): top panel in main
# Model used in the inference: homogeneous-response without differential mutant fitness (HOM0), 
# for purpose of comparison with previous results, as this model would have been used in previous studies
p_M_HOM0 <- ggplot(data = df, aes(x=ID, y=M_HOM0.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic), shape = 17) +
  geom_errorbar(aes(ymin=M_HOM0.2, ymax=M_HOM0.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + #guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle("Estimated using model HOM0")
p_M_HOM0
# Print experiments for which an increase in mutation rate is estimated (model = homogeneous-response without differential mutant fitness)
print(length(subset(df, M_HOM0.1 > 1)$ID))

# Fig1 version (ii): panel in SI
# Model used in the inference: the homogeneous-response/null model with the lowest AIC (N0, N1, N2, HOM0, HOM1 or HOM2)
p_M_AIC <- ggplot(data = df, aes(x=ID, y=M_AIC.1, group=antibiotic, shape = null_hom_by_AIC)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +  scale_shape_manual(values = c(2,5,0,17,18,15), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle("Estimates from the homogeneous model with lowest AIC")
p_M_AIC
# Print experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC)
print(length(subset(df, M_AIC.1 > 1)$ID))

# Fig1 version (iii): panel in SI
# Model used in the inference: the homogeneous-response/null model with the lowest AIC_corr
# AIC corrected for the sample size n: AIC_corr := 2p n/(n-p-1) - 2ln(L)
p_M_AIC_corr <- ggplot(data = df, aes(x=ID, y=M_AIC_corr.1, group=antibiotic, shape = null_hom_by_AIC_corr)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle(TeX("Estimates from the homogeneous model with lowest $AIC_c$"))
p_M_AIC_corr
# Print experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC_corr)
print(length(subset(df, M_AIC_corr.1 > 1)$ID))

# Null vs. HOM SI Figure version (i)
# Delta AIC_corr between the best model with M>1 (HOM0, HOM1 or HOM2) and the best null model (N0, N1, N2)
p_Delta_AIC_corr_NH <- ggplot(data = df, aes(x=ID, y=Delta_AIC_corr_NH, group=antibiotic)) + geom_point(aes(color=antibiotic, shape=by_AIC_corr_NH), size = 2) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  scale_y_continuous(limits = c(-10, 27)) + geom_point(data = subset(df, Delta_AIC_corr_NH > 27), aes(x=ID, y=Inf, color=antibiotic, shape=by_AIC_corr_NH), size=2) +
  geom_hline(yintercept = 2, linetype = "dashed") + geom_hline(yintercept = -2, linetype = "dashed") + 
  ylab(TeX("Difference in $AIC_c$")) + xlab("Experiment ID") + scale_shape_manual(values = c(17,18,15,2,5,0), name = "Overall selected") 
p_Delta_AIC_corr_NH
# For some experiments, Delta_AIC_corr is extremely large and therefore not shown in this plot
print(as.character(subset(df, Delta_AIC_corr_NH > 27)$ID))  # Experiment ID
print(subset(df, Delta_AIC_corr_NH > 27)$Delta_AIC_corr_NH) # Delta_AIC_corr


# Main/SI Figure: Experimental design impacting the results

# Do the width of the confidence intervals of the estimated change in population-wide mutation rate M depend on the experimental desgin, 
# i.e. plated fraction and number of parallel cultures
# Print: how many experiments use a plated fraction smaller than 1
print(length(subset(df, plated_fraction < 1)$ID))
df$width_CI_HOM0 <- (df$M_HOM0.3-df$M_HOM0.2)/df$M_HOM0.1
# Fitting a linear model with the log-width of the confidence interval of M estimated by model HOM0 as predicted variable
# and the log-plated fraction and log-number of parallel cultures as fixed effects
lm_HOM0 <- lm(data = df, log10(width_CI_HOM0) ~ log10(n_cultures) + log10(plated_fraction))
summary(lm_HOM0)
p_CI_HOM0 <- ggplot(data = df, aes(x=(plated_fraction), y=(width_CI_HOM0))) + #+log10(n_cultures)
  geom_point(aes(color=log10(n_cultures))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans="log10") +
  labs(x="Plated fraction", y="Normalised width of 95% CI around MLE estimate of M", color=TeX("$log_{10}(c_s)$")) +
  ggtitle("M estimated using model HOM0")
p_CI_HOM0
# Correlation between the plated fraction and the number of parallel cultures
cor.test(df$plated_fraction, df$n_cultures, method = "kendall")

# Same as above but with the CIs of M estimated by the homogeneous model with the lowest AIC_corr
df$width_CI_AIC_corr <- (df$M_AIC_corr.3-df$M_AIC_corr.2)/df$M_AIC_corr.1
lm_AIC_corr <- lm(data = subset(df, width_CI_AIC_corr > 0), log10(width_CI_AIC_corr) ~ log10(n_cultures) + log10(plated_fraction))
summary(lm_AIC_corr)
p_CI_AIC_corr <- ggplot(data = subset(df, width_CI_AIC_corr > 0), aes(x=(n_cultures), y=(width_CI_AIC_corr))) + #+log10(n_cultures)
  geom_point(aes(color=log10(plated_fraction))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans="log10") +
  labs(x="Number of parallel cultures", y="Normalised width of 95% CI around MLE estimate of M", color=TeX("$log_{10}(E)$")) +
  ggtitle(TeX("M estimates from the homogeneous model with lowest $AIC_c$"))
p_CI_AIC_corr
# Correlation between the plated fraction and the number of parallel cultures
cor.test(subset(df, width_CI_AIC_corr > 0)$plated_fraction, subset(df, width_CI_AIC_corr > 0)$n_cultures, method = "kendall")


# Goodness-of-fit test: could the data observed in the S condition have been generated under the parameters inferred from the UT condition?
# Print: experiments for which the S condition is significant different from the UT condition (criterion: p<0.05 in GoF)
print(length(subset(df, p_value_test_min < 0.05)$ID))


# Restriction to experiments for which we S is significantly different from UT
# (setting the estimated change in mutation rate to M=1, and the estimated supply-ration to S=0, if not)
df_strict <- df
print(as.character(subset(df_strict, p_value_test_min >= 0.05)$ID))
# Set M estimated by HOM0 to 1 if S/UT not significantly different
df_strict$M_HOM0.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_HOM0.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_HOM0.3[df_strict$p_value_test_min >= 0.05] <- 1.
# Creating new column to indicate whether S/UT significantly different (and we continue using HOM0) or not
df_strict$H0 <- "HOM0"
df_strict$H0[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
# Setting the Delta_AIC (and Delta_AIC_corr) between null and homogeneous model with M>1 to -Inf if S/UT is not significantly different
df_strict$Delta_AIC_NH[df_strict$p_value_test_min >= 0.05] <- -Inf
df_strict$by_AIC_NH[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
df_strict$Delta_AIC_corr_NH[df_strict$p_value_test_min >= 0.05] <- -Inf
df_strict$by_AIC_corr_NH[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
# Setting M=1 for the model with the lowest AIC (and AIC_corr) if S/UT not significantly
df_strict$M_AIC.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.3[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.2[df_strict$p_value_test_min >= 0.05] <- 1.
levels(df_strict$null_hom_by_AIC) <- c(levels(df_strict$null_hom_by_AIC), "S/UT not s. different")
df_strict$null_hom_by_AIC[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
df_strict$M_AIC_corr.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.3[df_strict$p_value_test_min >= 0.05] <- 1.
levels(df_strict$null_hom_by_AIC_corr) <- c(levels(df_strict$null_hom_by_AIC_corr), "S/UT not s. different")
df_strict$null_hom_by_AIC_corr[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
# And setting S=0 for the model with the lowest AIC (and AIC_corr) if S/UT not significantly
df_strict$S_AIC.1[df_strict$p_value_test_min >= 0.05] <- 0.
df_strict$S_AIC.3[df_strict$p_value_test_min >= 0.05] <- 0.
df_strict$S_AIC.2[df_strict$p_value_test_min >= 0.05] <- 0.
levels(df_strict$null_het_by_AIC) <- c(levels(df_strict$null_het_by_AIC), "S/UT not s. different")
df_strict$null_het_by_AIC[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"
df_strict$S_AIC_corr.1[df_strict$p_value_test_min >= 0.05] <- 0.
df_strict$S_AIC_corr.2[df_strict$p_value_test_min >= 0.05] <- 0.
df_strict$S_AIC_corr.3[df_strict$p_value_test_min >= 0.05] <- 0.
levels(df_strict$null_het_by_AIC_corr) <- c(levels(df_strict$null_het_by_AIC_corr), "S/UT not s. different")
df_strict$null_het_by_AIC_corr[df_strict$p_value_test_min >= 0.05] <- "S/UT not s. different"

# Evaluate the model fits using goodness-of-fit tests
# Print: experiments for which the homogeneous model with the lowest AIC is a poor fit to the data (none as of 21/07/2025)
print(as.character(subset(df_strict, p_value_UT_hom_AIC < 0.05)$ID)) # AIC in the UT condition
print(as.character(subset(df_strict, p_value_S_hom_AIC < 0.05)$ID))  # AIC in the S condition
print(as.character(subset(df_strict, p_value_UT_hom_AIC_corr < 0.05)$ID)) # AIC_corr in the UT condition
print(as.character(subset(df_strict, p_value_S_hom_AIC_corr < 0.05)$ID))  # AIC_corr in the S condition

# Same as Fig1 version (i), but marking experiments with S/UT not significantly different
p_M_HOM0_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_HOM0.1, group=antibiotic, shape = factor(H0))) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_HOM0.2, ymax=M_HOM0.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(17,8), name = "GoF test resulted in") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle("Estimated using model HOM0, if S significantly different from UT condition")
p_M_HOM0_GoF

# Same as Fig1 version (ii) but marking experiments with S/UT not significantly different
p_M_AIC_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_AIC.1, group=antibiotic, shape = null_hom_by_AIC)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15,8), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle("Estimates from the homogeneous model with lowest AIC, if S significantly different from UT condition")
p_M_AIC_GoF

# Same as Fig1 version (iii) but marking experiments with S/UT not significantly different, or poor model fits
p_M_AIC_corr_GoF <- ggplot(data = df_strict, aes(x=ID, y=M_AIC_corr.1, group=antibiotic, shape = null_hom_by_AIC_corr)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) + scale_shape_manual(values = c(2,5,0,17,18,15,8), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10", limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  ggtitle(TeX("Estimates from the homogeneous model with lowest $AIC_c$, if S significantly different from UT condition"))
p_M_AIC_corr_GoF

# Null vs. HOM SI Figure version (ii)
# Delta AIC_corr between the best model with M>1 (HOM0, HOM1 or HOM2) and the best null model (N0, N1, N2)
# Taking into account GoF test whether S/UT are significanlty different
p_Delta_AIC_corr_GoF_NH <- ggplot(data = df_strict, aes(x=ID, y=Delta_AIC_corr_NH, group=antibiotic)) + geom_point(aes(color=antibiotic, shape=by_AIC_corr_NH), size = 2) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_strict$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  scale_y_continuous(limits = c(-10, 27)) + geom_point(data = subset(df_strict, Delta_AIC_corr_NH > 27), aes(x=ID, y=Inf, color=antibiotic, shape=by_AIC_corr_NH), size=2) +
  geom_hline(yintercept = 2, linetype = "dashed") + geom_hline(yintercept = -2, linetype = "dashed") + 
  ylab(TeX("Difference in $AIC_c$")) + xlab("Experiment ID") + scale_shape_manual(values = c(17,18,15,2,5,0,8), name = "Overall selected") 
p_Delta_AIC_corr_GoF_NH
# Same as before, for some experiments Delta_AIC_corr is extremely large and therefore not shown in this plot
print(as.character(subset(df, Delta_AIC_corr_NH > 27)$ID))  # Experiment ID
print(subset(df, Delta_AIC_corr_NH > 27)$Delta_AIC_corr_NH) # Delta_AIC_corr

# Estimation of all experiments using heterogeneous-response models
# Evaluate the heterogeneous model fits using goodness-of-fit tests
# Print: E. coli experiments for which the heterogeneous model with the lowest AIC_corr is a poor fit to the data
print(as.character(subset(df_strict, p_value_UT_het_AIC_corr < 0.05)$ID)) # AIC_corr in the UT condition
print(as.character(subset(df_strict, p_value_S_het_AIC_corr < 0.05)$ID))  # AIC_corr in the S condition

# Generate a separate data frame to mark experiments for which the heterogeneous model fits under UT/S conditions are poor
df_strict_GoF_AIC_corr <- subset(df_strict, p_value_UT_het_AIC_corr < 0.05 | p_value_S_het_AIC_corr < 0.05)
levels(df_strict_GoF_AIC_corr$null_het_by_AIC_corr) <- c(levels(df_strict_GoF_AIC_corr$null_het_by_AIC_corr), "Poor fit UT cond", "Poor fit S cond")
df_strict_GoF_AIC_corr$null_het_by_AIC_corr[df_strict_GoF_AIC_corr$p_value_UT_het_AIC_corr < 0.05] <- "Poor fit UT cond"
df_strict_GoF_AIC_corr$null_het_by_AIC_corr[df_strict_GoF_AIC_corr$p_value_S_het_AIC_corr < 0.05] <- "Poor fit S cond"

print(as.character(subset(df_strict, M_HOM0.1 < 1)$ID))

# Equivalent to the figure above, but for the heterogeneous models
p_S_AIC_corr_GoF <- ggplot(data = df_strict, aes(x=ID, y=S_AIC_corr.1, group=antibiotic, shape = null_het_by_AIC_corr)) + 
  geom_point(aes(color=antibiotic)) + geom_errorbar(aes(ymin=S_AIC_corr.2, ymax=S_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1, linetype = "dashed") + scale_shape_manual(values = c(2,1,16,8,3,4), name = "Selected model") +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antimicrobial (abbr)") + 
  scale_y_continuous(trans="log10") + #, limits = c(1.4*10^-4, 7.2*10^5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Mutation-supply ratio") + xlab("Experiment ID") + 
  ggtitle(TeX("Estimates from the heterogeneous model with lowest $AIC_c$, if S significantly different from UT condition")) +
  geom_point(data = df_strict_GoF_AIC_corr, aes(x=ID, y=S_AIC_corr.1, shape=null_het_by_AIC_corr), size = 3)
p_S_AIC_corr_GoF

# Model selection between homogeneous and heterogeneous-response models
# Only plot experiments for which S/UT is significantly different
df_SIM_1 <- subset(df_strict, p_value_test_min < 0.05) 
# Exclude experiments for which the null model is selected both within homogeneous and heterogeneous models
df_SIM_1 <- subset(df_SIM_1, !(is.element(null_hom_by_AIC_corr, c("N0", "N1", "N2")) & null_het_by_AIC_corr == "N0"))

# Experiments where a differential mutant fitness is an alternative explanation to heterogeneity in mutation rates
print(as.character(subset(df_SIM_1, null_hom_by_AIC_corr == "N1")$ID)) # N1 is the homogeneous model with the lowest AIC_corr
print(as.character(subset(df_SIM_1, by_AIC_corr == "N1")$ID))          # N1 is selected overall
print(as.character(subset(df_SIM_1, null_hom_by_AIC_corr == "N2")$ID)) # N2 is the homogeneous model with the lowest AIC_corr
print(as.character(subset(df_SIM_1, by_AIC_corr == "N2")$ID))          # N2 is selected overall
print(as.character(subset(df_SIM_1, is.element(null_hom_by_AIC_corr, c("N1","N2")))$ID)) # N1 or N2 is the homogeneous model with the lowest AIC_corr
# Experiments where M>1 but no evidence of heterogeneity in mutation rates, i.e. N0 selected over HET0 and HET2
print(as.character(subset(df_SIM_1, null_het_by_AIC_corr == "N0")$ID))

# Reorder the IDs according to the different categories described above, for plotting
df_SIM_1$ID <- factor(df_SIM_1$ID, levels = as.character(c(subset(df_SIM_1, null_het_by_AIC_corr=="N0")$ID, subset(df_SIM_1, null_het_by_AIC_corr!="N0" & !is.element(null_hom_by_AIC_corr, c("N1", "N2")))$ID, subset(df_SIM_1, is.element(null_hom_by_AIC_corr, c("N1", "N2")))$ID)))

# Difference in AIC_corr between selected homogeneous and heterogeneous-response models including null models with differential mutant fitness N1 and N2
p_Delta_AIC_corr <- ggplot(data = df_SIM_1, aes(x=ID, y=Delta_AIC_corr, group=antibiotic)) + geom_point(aes(color=antibiotic, shape=by_AIC_corr), size = 2) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM_1$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  geom_hline(yintercept = 2, linetype = "dashed") + geom_hline(yintercept = -2, linetype = "dashed") + 
  ylab(TeX("Difference in $AIC_c$")) + xlab("Experiment ID") + scale_shape_manual(values = c(1,17,18,15,5,0,8,3,4), name = "Overall selected") 
p_Delta_AIC_corr

# Exclude the experiments where there no evidence for an increase M>1 overall
# either because the overall selected model is a null model
print(as.character(subset(df_SIM_1, is.element(by_AIC_corr, c("N0", "N1", "N2")))$ID))
df_SIM_1 <- subset(df_SIM_1, !is.element(by_AIC_corr, c("N0", "N1", "N2")))
# or a homogeneous model is selected overall which estimates M<1
print(as.character(subset(df_SIM_1, (is.element(by_AIC_corr, c("HOM0", "HOM1", "HOM2")) & M_AIC_corr.1 < 1) & !is.element(by_AIC_corr, c("HET0", "HET2")))$ID))
df_SIM_1 <- subset(df_SIM_1, (is.element(by_AIC_corr, c("HOM0", "HOM1", "HOM2")) & M_AIC_corr.1 > 1) | is.element(by_AIC_corr, c("HET0", "HET2")))

# Define 'stress-induced mutagenesis is detected' for homogeneous models only (from least to most strict)
# 1. The homogeneous model with the lowest AIC estimates an increase in population-wide mutation rate M>1 (the MLE, not the bounds of the CI)
df$SIM_AIC <- df$M_AIC.1 > 1
print(length(subset(df, SIM_AIC == TRUE)$ID))
# 2. The homogeneous model with the lowest AIC_corr estimates an increase M>1
df$SIM_AIC_corr <- df$M_AIC_corr.1 > 1
print(length(subset(df, SIM_AIC_corr == TRUE)$ID))
# 3. The homogeneous model with the lowest AIC estimates an increase M>1, and the S condition is significantly different from the UT condition
df_strict$SIM_AIC <- df_strict$M_AIC.1 > 1
print(length(subset(df_strict, SIM_AIC == TRUE)$ID))
# 4. The homogeneous model with the lowest AIC_corr estimates an increase M>1, and S/UT significantly different
df_strict$SIM_AIC_corr <- df_strict$M_AIC_corr.1 > 1
print(length(subset(df_strict, SIM_AIC_corr == TRUE)$ID))

# Using definition 4, but with this definition of SIM we exclude cases in which HET0 is selected over N1/N2 or HOM with M<1
df_SIM_2 <- subset(df_strict, SIM_AIC_corr == TRUE)
print(length(df_SIM_2$ID)) # Only experiments with M>1
print(length(df_SIM_1$ID)) # Including experiments for which HET0 is selected over models with M<=1
print(as.character(setdiff(df_SIM_1$ID, df_SIM_2$ID)))
print(as.character(setdiff(subset(df_SIM_1, null_hom_by_AIC_corr == "N0")$ID, df_SIM_2$ID))) # N0 is selected over homogeneous models
print(as.character(setdiff(subset(df_SIM_1, null_hom_by_AIC_corr == "N1")$ID, df_SIM_2$ID))) # N1 is selected over homogeneous models
print(as.character(setdiff(subset(df_SIM_1, null_hom_by_AIC_corr == "N2")$ID, df_SIM_2$ID))) # N2 is selected over homogeneous models

# Comparing to for how many experiments SIM would have been detected solely on the basis of estimating an increase M>1
print(length(subset(df, M_HOM0.1 > 1)$ID))

# Estimated M for experiments with detected SIM in any form (homogeneous or heterogeneous)
# Estimates for which a heterogeneous is selected over a homogeneous response overall are shown in dashed lines
df_SIM_1$homhet <- character(length(df_SIM_1$by_AIC_corr))
df_SIM_1$homhet[is.element(df_SIM_1$by_AIC_corr, c("N1", "N2", "HOM0", "HOM1", "HOM2"))] <- "Homogeneous"
df_SIM_1$homhet[is.element(df_SIM_1$by_AIC_corr, c("N0", "HET0", "HET2"))] <- "Heterogeneous"
df_SIM_1$homhet_sel <- df_SIM_1$null_hom_by_AIC_corr == df_SIM_1$by_AIC_corr
p_M_SIM <- ggplot(data = df_SIM_1, aes(x=ID, y=M_AIC_corr.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic, shape = null_hom_by_AIC_corr), size=2) + scale_shape_manual(values = c(2,5,17,18,15), name = "Selected model") +
  geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic, linetype = homhet_sel)) + geom_hline(yintercept = 1) +
  scale_colour_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM_1$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Overall selected") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate") + xlab("Experiment ID")
p_M_SIM

# Estimated M for experiments with detected SIM in any form (homogeneous or heterogeneous)

# Again, generate a separate data frame to mark experiments for which the heterogeneous model fits under UT/S conditions are poor
df_SIM_GoF_AIC_corr <- subset(df_SIM_1, p_value_UT_het_AIC_corr < 0.05 | p_value_S_het_AIC_corr < 0.05)
levels(df_SIM_GoF_AIC_corr$null_het_by_AIC_corr) <- c(levels(df_SIM_GoF_AIC_corr$null_het_by_AIC_corr), "Poor fit UT cond", "Poor fit S cond")
df_SIM_GoF_AIC_corr$null_het_by_AIC_corr[df_SIM_GoF_AIC_corr$p_value_UT_het_AIC_corr < 0.05] <- "Poor fit UT cond"
df_SIM_GoF_AIC_corr$null_het_by_AIC_corr[df_SIM_GoF_AIC_corr$p_value_S_het_AIC_corr < 0.05] <- "Poor fit S cond"

# Estimates for which a homogeneous is selected over a heterogeneous response overall are shown in dashed lines
p_S_SIM <- ggplot(data = df_SIM_1, aes(x=ID, y=S_AIC_corr.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic, shape = null_het_by_AIC_corr), size=2) + scale_shape_manual(values = c(2,1,16,3,4), name = "Selected model") +
  geom_errorbar(aes(ymin=S_AIC_corr.2, ymax=S_AIC_corr.3, color=antibiotic, linetype = !homhet_sel)) +
  scale_linetype_manual(values = c("dashed", "solid"), name = "Overall selected") + geom_hline(yintercept = 1, linetype = "dashed") +
  scale_colour_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM_1$antibiotic)))$color, name = "Antimicrobial (abbr)") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Mutation-supply ratio") + xlab("Experiment ID") + geom_point(data = df_SIM_GoF_AIC_corr, aes(x=ID, y=S_AIC_corr.1, group=antibiotic, shape=null_het_by_AIC_corr), size = 2)
p_S_SIM

# We restrict the analysis about the target to experiments using E. coli MG1655 and TD2158 (no mutant strains)
df_Ecoli <- subset(df, species == "E. coli" & is.element(strain, c("MG1655", "TD2158")))
df_Ecoli_strict <- subset(df_strict, species == "E. coli" & is.element(strain, c("MG1655", "TD2158")))

# How many E. coli experiments use DNA/DNA gyrase-targeting antibiotics?
length(subset(df_Ecoli, target_group == "DNA/DNA gyrase")$ID)
# And how many use ribosome-targeting ones?
length(subset(df_Ecoli, target_group == "Ribosome")$ID)
# And how many use others?
length(subset(df_Ecoli, target_group == "Other")$ID)

# Testing for normality of the estimated increase M
shapiro.test(df_Ecoli$M_HOM0.1)                                           # All E. coli experiments
shapiro.test(subset(df_Ecoli, target_group == "DNA/DNA gyrase")$M_HOM0.1) # E. coli experiments using a DNA/DNA gyrase-targeting antibiotic
shapiro.test(subset(df_Ecoli, target_group == "Ribosome")$M_HOM0.1)       # E. coli experiments using a ribosome-targeting antibiotic
shapiro.test(subset(df_Ecoli, target_group == "Other")$M_HOM0.1)          # E. coli experiments using another antibiotic

# Kruskal-Wallis test -> Does the estimated M depend on the antibiotic target?
kruskal.test(M_HOM0.1 ~ target, data = df_Ecoli)
# More specifically, is M significantly different for DNA/DNA gyrase vs. ribosome-targeting antibiotics or others?
kruskal.test(M_HOM0.1 ~ target_group, data = df_Ecoli)
median(subset(df_Ecoli, target_group == "DNA/DNA gyrase")$M_HOM0.1)
median(subset(df_Ecoli, target_group == "Ribosome")$M_HOM0.1)
median(subset(df_Ecoli, target_group == "Other")$M_HOM0.1)
# Follow-up pairwise test: Dunn test with Holm p-value adjustment 
dunnTest(M_HOM0.1 ~ as.factor(target_group), data = df_Ecoli, method = "holm")
#pairwise.wilcox.test( df_Ecoli$M_HOM0.1, df_Ecoli$target_group, p.adjust.method = "bonferroni")

# We also plot experiments for which SIM is detected here (defined as model with lowest AIC_corr estimates M>1 and S/UT significantly different)
# For how many experiments do we detect SIM, depending on the target?
sum(subset(df_Ecoli_strict, target_group == "DNA/DNA gyrase")$SIM_AIC_corr) # DNA/DNA gyrase-targeting 
sum(subset(df_Ecoli_strict, target_group == "Ribosome")$SIM_AIC_corr)       # Ribosome-targeting
sum(subset(df_Ecoli_strict, target_group == "Other")$SIM_AIC_corr)          # Others
# For this plot, I want the y value to be the estimate using HOM0 without taking the GoF test into account (to avoid repeated values y=1)
# and the colour according to our definition of SIM which takes the GoF test into account
df_Ecoli_plot <- df_Ecoli
df_Ecoli_plot$SIM_AIC_corr <- df_Ecoli_strict$SIM_AIC_corr
p_M_DNA <- ggplot(data = df_Ecoli, aes(x=target_group, y=M_HOM0.1)) + geom_boxplot(aes(fill=target_group), show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(data = df_Ecoli_plot, aes(color = SIM_AIC_corr), width = 0.25, alpha = 0.8) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey"), name = "SIM") +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  scale_fill_manual(values = c("DNA/DNA gyrase" = all_colours[round(n_DNA)], "Other" = all_colours[round(n_DNA+n_Other)], "Ribosome" = all_colours[round(n_DNA+n_Other+n_Ribosome)])) + 
  theme(plot.margin = margin(0.5,0.5,2.5,0.5, "cm")) + 
  ylab("Fold change in population-wide mutation rate estimated using model HOM0") + xlab("Antimicrobial target") +
  stat_compare_means(method = "kruskal.test", label.y = 300) + theme(legend.position = "right")
p_M_DNA
print(chisq.test(df_Ecoli_strict$target_group, df_Ecoli_strict$SIM_AIC_corr))


# What else does detection of SIM depend on? Target/experimental design?
# For our linear mixed model, we again define SIM as 'homogeneous model with the lowest AIC_corr estimates M>1 and S/UT significantly different'
# and restrict the analysis to experiments using E. coli wt strains
df_glm <- df_Ecoli_strict
#df_glm <- subset(df_Ecoli_strict, target_group != "Ribosome") # -> no difference between 'Other' and 'DNA/DNA gyrase-targeting
df_glm$target_group <- relevel(factor(df_glm$target_group), ref = "Ribosome")
# Fitting a linear mixed model with detection of SIM as a predicted variable
# and the target group and the log-plated fraction as fixed effects
glm <- glm(SIM_AIC_corr ~ target_group + log10(n_cultures) + log10(plated_fraction), data = df_glm)
summary(glm)

# Previously, I tried fitting a generalised linear mixed model with 
# the antibiotic concentration in units of the MIC as another fixed effect
# and the ID as a random effect (here, baseline_ID as a random effect does not account for any variance)
# Problems: for the E. coli experiments with known MIC alone complete separate occurred because no SIM is detected for any ribosome-targeting antibiotic
df_glmm <- df_Ecoli_strict
df_glmm$target_group <- relevel(factor(df_glmm$target_group), ref = "Ribosome")
df_glmm <- subset(df_glmm, !is.na(of_MIC))
length(df_Ecoli_strict$ID) # All E. coli experiments
length(df_glmm$ID)         # Experiments for which the MIC was measured
glmm <- glmer(SIM_AIC_corr ~ of_MIC + target_group + log10(n_cultures) + log10(plated_fraction) + (1|ID), data = df_glmm, family = binomial)
summary(glmm)
# Problems did not occur when using experiments with all species
# Note that using the baseline_ID as a random effect should deal with difference between species and strains 
df_glmm <- df_strict
df_glmm$target_group <- relevel(factor(df_glmm$target_group), ref = "Ribosome")
df_glmm <- subset(df_glmm, !is.na(of_MIC))
glmm <- glmer(SIM_AIC_corr ~ of_MIC + target_group + log10(n_cultures) + log10(plated_fraction) + (1|baseline_ID), data = df_glmm, family = binomial)
summary(glmm)

# Overall model selection (criterion = lowest AIC_corr) for all E. coli wt experiments
# 3 categories: 'null' (models N0, N1 or N2), 'homogeneous' (models HOM0, HOM1 or HOM2) and 'heterogeneous' (models HET0 or HET2)
# Exclude experiments for which a homogeneous model with estimated M<1 is selected
df_Ecoli_sel <- subset(df_Ecoli_strict, !(is.element(by_AIC_corr, c("HOM0", "HOM1", "HOM2")) & M_AIC_corr.1 < 1))
df_sel <- data.frame(antibiotic=rep(unique(df_Ecoli_sel$antibiotic), each=3))
df_sel$target <- plyr::mapvalues(df_sel$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target))
df_sel$target_group <- plyr::mapvalues(df_sel$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target_group))
df_sel$nullhomhet <- factor(rep(c("Null","Homogeneous","Heterogeneous"), length(unique(df_Ecoli_sel$antibiotic))), levels = c("Null", "Homogeneous", "Heterogeneous"), ordered = TRUE)
v <- numeric(length(df_sel$antibiotic))
for (i in 1:length(unique(df_Ecoli_sel$antibiotic))) {
  v[3*i-2] <- sum(is.element(subset(df_Ecoli_sel, antibiotic == df_sel$antibiotic[3*i])$by_AIC_corr, c("N0","N1","N2"))) 
  v[3*i-1] <- sum(is.element(subset(df_Ecoli_sel, antibiotic == df_sel$antibiotic[3*i])$by_AIC_corr, c("HOM0","HOM1","HOM2"))) 
  v[3*i] <- sum(is.element(subset(df_Ecoli_sel, antibiotic == df_sel$antibiotic[3*i])$by_AIC_corr, c("HET0","HET2"))) 
}
df_sel$prevalence <- v
# Grouping by antibiotic target group ('DNA/DNA gyrase', 'Other' or 'Ribosome')
p_sel_t <- ggplot(data = df_sel, aes(x=nullhomhet, y=prevalence, fill=target_group)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("DNA/DNA gyrase" = all_colours[round(n_DNA)], "Other" = all_colours[round(n_DNA+n_Other)], "Ribosome" = all_colours[round(n_DNA+n_Other+n_Ribosome)]), name = "Target") +
  xlab("Overall selected") + ggtitle(TeX("Selection criterion: lowest $AIC_c$")) + ylab("Number of experiments") 
p_sel_t
# How many experiments are in each of the categories? ('DNA/DNA gyrase', 'Other' or 'Ribosome')
# Null model selected
print(c(sum(subset(df_sel, target_group=="DNA/DNA gyrase" & nullhomhet=="Null")$prevalence), sum(subset(df_sel, target_group=="Other" & nullhomhet=="Null")$prevalence), sum(subset(df_sel, target_group=="Ribosome" & nullhomhet=="Null")$prevalence))) 
# Homogeneous model selected
print(c(sum(subset(df_sel, target_group=="DNA/DNA gyrase" & nullhomhet=="Homogeneous")$prevalence), sum(subset(df_sel, target_group=="Other" & nullhomhet=="Homogeneous")$prevalence), sum(subset(df_sel, target_group=="Ribosome" & nullhomhet=="Homogeneous")$prevalence))) 
# Heterogeneous model selected
print(c(sum(subset(df_sel, target_group=="DNA/DNA gyrase" & nullhomhet=="Heterogeneous")$prevalence), sum(subset(df_sel, target_group=="Other" & nullhomhet=="Heterogeneous")$prevalence), sum(subset(df_sel, target_group=="Ribosome" & nullhomhet=="Heterogeneous")$prevalence))) 


# Case studies
# Experiments for which the selected homogeneous/heterogeneous models have comparable AIC_corrs
print(as.character(subset(df_Ecoli_sel, Delta_AIC_corr > 0)$ID))
# Experiments for which the homogeneous model with unconstrained mutant fitness is selected
print(as.character(subset(df_Ecoli_sel, null_hom_by_AIC_corr == "HOM2")$ID))


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

df_Nor <- subset(est_paras, ID == "Frenoy_Nor")
# UT condition: calculate the mutant count distributions under HOM2, HOM1, HET2 and HET2 with mutant fitness as estimated by HOM2
p_mc_HOM2_UT <- pMudi(max(mc_UT), Nf_UT, paras_HOM2$MLE[1], plateff=eff_UT, fit_m=paras_HOM2$MLE[2]) * length(mc_UT)
p_mc_HOM1_UT <- pMudi(max(mc_UT), Nf_UT, paras_HOM1$MLE[1], plateff=eff_UT, fit_m=paras_HOM1$MLE[2]) * length(mc_UT)
p_mc_HET2_UT <- pMudi(max(mc_UT), Nf_UT, paras_HET2$MLE[1], plateff=eff_UT) * length(mc_UT)
p_mc_HET2_fitm_UT <- pMudi(max(mc_UT), Nf_UT, paras_HET2_fitm$MLE[1], plateff=eff_UT, fit_m=paras_HET2_fitm$MLE[2]) * length(mc_UT)
# Estimated mutation rate and mutant fitness under HOM2, and estimated mutation rate and mutant fitness under HOM1
print(c(paras_HOM2$MLE[1], paras_HOM2$MLE[2], paras_HOM1$MLE[1], paras_HOM1$MLE[2]))
# Estimated mutation rate off-cells and fixed mutant fitness under HET2, and estimated mutation rate off-cells and fixed mutant fitness under HET2 with mutant fitness from HOM2
print(paras_HET2$MLE[1], paras_HET2$MLE[2], paras_HET2_fitm$MLE[1], paras_HET2_fitm$MLE[2])
p_mc_UT <- ggplot() + geom_histogram(aes(mc_UT, y=after_stat(density) * length(mc_UT)), fill="#1C458A", bins = 50) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HOM2_UT), color="#F0E442", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HOM1_UT), color="#E69F00", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HET2_UT), color="#CC79A7", linewidth = 1., linetype = "dashed") + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_HET2_fitm_UT), color="#7570B3", linewidth = 1., linetype = "dashed") + 
  ggtitle("Untreated experiments") + xlab("Number of colonies") + ylab("Number of plates") + ylim(-0.1, 7.1)
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
  geom_line(aes(x=0:max(mc_S), y=p_mc_HET2_fitm_S), color="#7570B3", linewidth = 1., linetype = "dashed") +
  xlab("Number of colonies") + ylab("Number of plates") + ggtitle("Experiments with norfloxacin") + ylim(-0.5,34)
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
mc_data <- read.table(paste0("experimental_data/raw_counts/Baharoglu_Cip.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])       # Mutant counts
Nf_S <- mean(read_counts(mc_data[3,])) # Average final population size
eff_S <- as.numeric(mc_data[4,1])      # Plating efficiency

# Estimation under model HOM0 (by setting fit_m = c(1., 1.))
res <- estimu(mc_UT, Nf_UT, mc_S, Nf_S, plateff = c(eff_UT, eff_S), fit_m = c(1., 1.), mod = "hom")
# Data frame with the estimated parameters
paras_HOM0 <- res[[1]]
print(paras_HOM0)
# Data frame with the results relevant for model selection
msel_HOM0 <- res[[2]]
print(msel_HOM0)

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
paras_HOM0$MLE[6]     # Homogeneous (without differential mutant fitness)




# Further analysis

# Does it depend on the experimental design impact the outcome of the GoF to test whether S/UT are significantly different?
df_GoF <- subset(df, p_value_test_min > 0)
# Fitting a linear model with the p-value as predicted variable, and log-plated fraction and log-number of cultures as predictors
lm_p <- lm(p_value_test_min ~ log10(plated_fraction) + log10(n_cultures), data = df_GoF)
summary(lm_p)
p_GoF <- ggplot(data = df_GoF, aes(x=(plated_fraction), y=(p_value_test_min))) + #+log10(n_cultures)
  geom_point(aes(color=log10(n_cultures))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + #scale_y_continuous(trans = "log10") +
  labs(x="Plated fraction", y="p value GoF test (S generated under UT parameters)", color=TeX("$log_{10}(c_s)$"))
p_GoF

# Are p-values of the GoF tests different for the UT and S conditions?
# My intuition here is that in the S condition it should be more likely that something completely different is going on, 
# which cannot be explained well by any of the models, not even the homogeneous model with the lowest AIC_corr
df_p_value_UT_S <- data.frame(ID = rep(df$ID, 2))
df_p_value_UT_S$condition <- rep(c("UT", "S"), each = length(df$ID))
df_p_value_UT_S$p_values <- c(df$p_value_UT_hom_AIC_corr, df$p_value_S_hom_AIC_corr)
wilcox.test(p_values ~ condition, data = df_p_value_UT_S)
median(subset(df_p_value_UT_S, condition == "UT")$p_values)
median(subset(df_p_value_UT_S, condition == "S")$p_values)
p_M_DNA <- ggplot(data = df_p_value_UT_S, aes(x=condition, y=p_values)) + 
  geom_boxplot(aes(fill=condition), show.legend = FALSE, outlier.shape = NA) + 
  ylab("p-value goodness of fit test") + xlab("Treatment condition") + ggtitle(TeX("Hom model with lowest $AIC_c$")) + 
  stat_compare_means(label.y = 1) + theme(legend.position = "right")
p_M_DNA

# I also tried fitting a linear mixed model with the p-value of the S/UT GoF as a predicted variable 
# and the antibiotic dose, the target group and the log-plated fraction as fixed effects
# and the baseline ID as random effect -> doesn't account for any variance when restricting to E. coli wt experiments!
lm_p <- lm(p_value_test_min ~ of_MIC + target_group + log10(plated_fraction) + log10(n_cultures), data = df_glmm)
summary(lm_p)
lm_p_wo_MIC <- lm(p_value_test_min ~ target_group + log10(plated_fraction) + log10(n_cultures), data = df_glmm)
summary(lm_p_wo_MIC)

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