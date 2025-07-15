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

# Add target to meta data, pool replicates and add levels to data frames
meta_data$target <- mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$target_group)
meta_data <- subset(meta_data, replicate == 0)
df <- merge(meta_data, sum_data, by = "ID")
df <- merge(df, est_sum, by = "ID")
df$antibiotic <- factor(df$antibiotic, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
df$target <- factor(df$target, levels = unique(antibiotic_classes$target_group), ordered = TRUE)
antibiotic_classes$antibiotic_abbr <- factor(antibiotic_classes$antibiotic_abbr, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
antibiotic_classes$target_group <- factor(antibiotic_classes$target_group, levels = unique(antibiotic_classes$target_group), ordered = TRUE)

# Color-coding antibiotics used in the studies (sorted by target group)
antibiotic_classes$color <- turbo(length(antibiotic_classes$antibiotic_abbr))
v <- numeric(length(antibiotic_classes$antibiotic_abbr))
for (i in 1:length(v)) {
  v[i] <- length(subset(meta_data, antibiotic == antibiotic_classes$antibiotic_abbr[i])$antibiotic)
}
antibiotic_classes$prevalence <- v
p_antibiotic <- ggplot(data = antibiotic_classes, aes(x=target_group, y=prevalence, fill=antibiotic_abbr)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = antibiotic_classes$color, name = "Antibiotic (abbr)") + 
  xlab("Grouped by target") + ylab("Number of experiments") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.8))
p_antibiotic

# Sort by bacterial species
df <- arrange(df, species, target)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)

# Goodness of fit test for the untreated condition
print(as.character(subset(df, p_value_UT_max < 0.05)$ID))
df <- subset(df, p_value_UT_max >= 0.05)

# Estimated increase in population-wide mutation rate by antibiotic
# Under homogeneous-response model without differential mutant fitness used in the inference (standard model for UT and S separately), 
# for purpose of comparison
p_M_antibiotic <- ggplot(data = df, aes(x=ID, y=M_wo_fitm.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_wo_fitm.2, ymax=M_wo_fitm.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

# Under homogeneous-response model with the lowest AIC
p_M_antibiotic <- ggplot(data = df, aes(x=ID, y=M_AIC.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

# Experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC)
print(length(subset(df, M_AIC.1 > 1)$ID))
# Experiments for which we reject hypothesis that S data was generated under UT parameters (criterion: p<0.05 in GoF)
print(length(subset(df, p_value_test_min < 0.05)$ID))

setdiff(subset(df, M_AIC.1 != 1)$ID, subset(df, p_value_test_min < 0.05)$ID)

# p value of GoF depending on experimental design
df_GoF <- df #subset(df, p_value_test_min > 0)
lm_p <- lm(p_value_test_min ~ log10(plated_fraction) + log10(n_cultures), data = df_GoF)
summary(lm_p)
p_GoF <- ggplot(data = df_GoF, aes(x=(n_cultures), y=(p_value_test_min))) + #+log10(n_cultures)
  geom_point(aes(color=log10(plated_fraction))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + #scale_y_continuous(trans = "log10") +
  labs(x="Number of parallel cultures", y="p value GoF test (S generated under UT parameters)", color=TeX("$log_{10}(E)$"))
p_GoF

# Restriction to experiments with rejected hypothesis 'S generated under UT parameters'
df_strict <- df
df_strict$M_AIC.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC.3[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.1[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.2[df_strict$p_value_test_min >= 0.05] <- 1.
df_strict$M_AIC_corr.3[df_strict$p_value_test_min >= 0.05] <- 1.
# 'Strict' model selection criterion: lowest AIC + GoF
df_strict$SIM <- df_strict$M_AIC.1 > 1
# Alternatively, using AIC_corr as the model selection criterion
df$SIM <- df$M_AIC_corr.1 > 1

p_M_antibiotic <- ggplot(data = df_strict, aes(x=ID, y=M_AIC.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

# Distribution of p-values of the goodness of fit test
est_paras_pooled <- subset(est_paras, is.element(ID, meta_data$ID))
crit <- "Hom with lowest AIC"
cond <- "UT condition"
p_values <- df$p_value_UT_hom_AIC
print(as.character(subset(df, p_value_S_het_AIC_corr < 0.05)$ID))
#p_values <- subset(est_paras_pooled, model == "het_div")$p_value_UT
res <- shapiro.test(p_values)
p_value <- res$p.value
hist(p_values, probability = TRUE, xlab = paste0("p-value ", cond), main = crit)
mtext(text = paste("Shapiro-Wilk p-value:", round(p_value, 3)), side = 3, line = 0.5, adj = 1)

# Are p-values different for the UT and S conditions?
df$SIM_lenient <- df$M_AIC.1 > 1
df_p_value_UT_S <- data.frame(ID = rep(df$ID, 2))
df_p_value_UT_S$condition <- rep(c("UT", "S"), each = length(df$ID))
df_p_value_UT_S$SIM_lenient <- rep(df$SIM_lenient, 2)
df_p_value_UT_S$SIM <- rep(df$SIM, 2)
crit <- "Hom model with lowest AIC_corr"
#mod <- "het_div"
df_p_value_UT_S$p_values <- c(df$p_value_UT_hom_AIC_corr, df$p_value_S_hom_AIC_corr)
#df_p_value_UT_S$p_values <- c(subset(est_paras_pooled, model == mod)$p_value_UT, subset(est_paras_pooled, model == mod)$p_value_S)
wilcox.test(p_values ~ condition, data = df_p_value_UT_S)
median(subset(df_p_value_UT_S, condition == "UT")$p_values)
median(subset(df_p_value_UT_S, condition == "S")$p_values)
p_M_DNA <- ggplot(data = df_p_value_UT_S, aes(x=condition, y=p_values)) + geom_boxplot(aes(fill=condition), show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(aes(color = SIM), width = 0.25, alpha = 0.8) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
  #coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  #scale_fill_manual(values = c("DNA/DNA gyrase" = "#4E6ADB", "Ribosome" = "#FE9B2D")) + 
  #theme(plot.margin = margin(0.5,0.5,2.5,0.5, "cm")) + 
  ylab("p-value goodness of fit test") + xlab("Treatment condition") + ggtitle(crit) + 
  stat_compare_means(label.y = 1) + theme(legend.position = "right")
p_M_DNA

# Under homogeneous-response model with the lowest AIC_corr
# AIC corrected for the sample size n: AIC_corr = 2p n/(n-p-1) - 2ln(L)
p_M_antibiotic <- ggplot(data = df, aes(x=ID, y=M_AIC_corr.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

# Experiments for which an increase in mutation rate is estimated (model selection criterion: lowest AIC_corr)
print(length(subset(df, M_AIC_corr.1 > 1)$ID))
# Experiments for which we reject hypothesis that S data was generated under UT parameters (criterion: p<0.05 in GoF)
print(length(subset(df, p_value_test_min < 0.05)$ID))

setdiff(subset(df, M_AIC_corr.1 != 1)$ID, subset(df, p_value_test_min < 0.05)$ID)
setdiff(subset(df, M_AIC.1 != 1)$ID, subset(df, M_AIC_corr.1 != 1)$ID)

p_M_antibiotic <- ggplot(data = df_strict, aes(x=ID, y=M_AIC_corr.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_AIC_corr.2, ymax=M_AIC_corr.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic


# Plating efficiency/number of parallel cultures and width of confidence intervals

print(length(subset(df, plated_fraction < 1)$ID))
df$width_CI <- (df$M_wo_fitm.3-df$M_wo_fitm.2)/df$M_wo_fitm.1
lm <- lm(log10(width_CI) ~ log10(plated_fraction) + log10(n_cultures), data = df)
summary(lm)
p_CI <- ggplot(data = df, aes(x=(plated_fraction), y=(width_CI))) + #+log10(n_cultures)
  geom_point(aes(color=log10(n_cultures))) + #geom_smooth(method = "lm") + 
  scale_x_continuous(trans = "log10") + scale_y_continuous(trans="log10") +
  labs(x="Plated fraction", y="Normalised width of 95% CI around MLE estimate", color=TeX("$log_{10}(c_s)$"))
p_CI

cor.test(df$plated_fraction, df$n_cultures, method = "kendall")

df_glmm <- df
# Model selection criterion for SIM: M > 1 for model with lowest AIC + GoF test 
df_glmm$SIM <- df_strict$M_AIC.1 > 1
df_glmm <- subset(df_glmm, !is.na(of_MIC))
df_glmm$group <- character(length(df_glmm$target))
df_glmm$group[is.element(df_glmm$target, c("DNA", "DNA gyrase"))] <- "DNA/DNA gyrase"
df_glmm$group[df_glmm$target == "Ribosome"] <- "Ribosome"
df_glmm$group[!is.element(df_glmm$target, c("DNA", "DNA gyrase", "Ribosome"))] <- "Other"
glmm <- glmer(SIM ~ of_MIC + group + log10(plated_fraction) + (1|baseline_ID), data = df_glmm, family=binomial)
summary(glmm)

lmer_p <- lmer(p_value_test_min ~ of_MIC + group + log10(plated_fraction) + log10(n_cultures) + (1|baseline_ID), data = df_glmm)
summary(lmer_p)

# Experiments with detected SIM
print(as.character(subset(df_strict, SIM == TRUE)$ID)) # AIC + GoF test
print(as.character(subset(df, SIM == TRUE)$ID))        # AIC_corr
# Our model selection procedure (AIC+GoF test, AIC_corr) compared to none
print(c(length(subset(df_strict, M_AIC.1 > 1)$ID), length(subset(df, M_AIC_corr.1 > 1)$ID), length(subset(df, M_wo_fitm.1 > 1)$ID)))

# Further analysis with experiments using E. coli MG1655 and TD2158 (no mutant strains)
df_Ecoli <- subset(subset(df, species == "E. coli"), is.element(strain, c("MG1655", "TD2158")))
df_Ecoli_strict <- subset(subset(df_strict, species == "E. coli"), is.element(strain, c("MG1655", "TD2158")))
print(as.character(subset(df_Ecoli_strict, M_AIC.1 > 1)$ID))
print(as.character(subset(df_Ecoli, M_AIC.1 > 1)$ID))

# Testing for normality -> not normal
shapiro.test(df_Ecoli$M_wo_fitm.1)
length(subset(df_Ecoli, is.element(target, c("DNA", "DNA gyrase")))$ID)
shapiro.test(subset(df_Ecoli, is.element(target, c("DNA", "DNA gyrase")))$M_wo_fitm.1)
length(subset(df_Ecoli, target=="Ribosome")$ID)
shapiro.test(subset(df_Ecoli, target == "Ribosome")$M_wo_fitm.1)
#hist(subset(df_Ecoli, target == "Ribosome")$M_wo_fitm.1, main="Ribosome-targeting antibiotics", xlab="Increase in population-wide mutation rate", probability=TRUE) +
#  abline(v = mean(subset(df_Ecoli, target == "Ribosome")$M_wo_fitm.1), col="red", lwd=2)

# Kruskal-Wallis test -> DNA/DNA-gyrase and ribosome binding significantly different
kruskal.test(M_wo_fitm.1 ~ target, data = df_Ecoli_strict)
df_KW <- subset(df_Ecoli_strict, is.element(target, c("DNA", "DNA gyrase", "Ribosome")))
df_KW$group <- character(length(df_KW$target))
df_KW$group[is.element(df_KW$target, c("DNA", "DNA gyrase"))] <- "DNA/DNA gyrase"
df_KW$group[df_KW$target == "Ribosome"] <- "Ribosome"
kruskal.test(df_KW$M_wo_fitm.1, df_KW$group)
wilcox.test(M_wo_fitm.1 ~ group, data = df_KW)
median(subset(df_KW, group == "DNA/DNA gyrase")$M_wo_fitm.1)
median(subset(df_KW, group == "Ribosome")$M_wo_fitm.1)
p_M_DNA <- ggplot(data = df_KW, aes(x=group, y=M_wo_fitm.1)) + geom_boxplot(aes(fill=group), show.legend = FALSE, outlier.shape = NA) + 
  geom_jitter(aes(color = SIM), width = 0.25, alpha = 0.8) + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkgrey")) +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  scale_fill_manual(values = c("DNA/DNA gyrase" = "#4E6ADB", "Ribosome" = "#FE9B2D")) + 
  theme(plot.margin = margin(0.5,0.5,2.5,0.5, "cm")) + 
  ylab("Estimated fold change in population-wide mutation rate") + xlab("Antimicrobial target") +
  stat_compare_means(label.y = 300) + theme(legend.position = "right")
p_M_DNA

print(chisq.test(df_KW$group, df_KW$SIM))

# Experiments for which SIM was detected
df_SIM <- subset(df_Ecoli_strict, SIM == TRUE)
df_SIM <- arrange(df_SIM, target)
df_SIM$ID <- factor(df_SIM$ID, levels = unique(df_SIM$ID), ordered = TRUE)
df_SIM$homhet <- character(length(df_SIM$by_AIC))
df_SIM$homhet[is.element(df_SIM$by_AIC, c("hom_wo_fitm", "hom_fitm", "hom_fitm_unconstr"))] <- "Homogeneous"
df_SIM$homhet[is.element(df_SIM$by_AIC, c("het_zero_div", "het_div"))] <- "Heterogeneous"
print(c(length(df_SIM$ID), length(subset(df_Ecoli, M_wo_fitm.1>1)$ID)))
print(c(length(subset(df_SIM, is.element(target, c("DNA", "DNA gyrase")))$ID),length(subset(df_SIM, target=="Ribosome")$ID)))
p_M_antibiotic <- ggplot(data = df_SIM, aes(x=ID, y=M_AIC.1, group=target)) + 
  geom_point(aes(color=target, shape = homhet), size=2) + scale_shape_manual(values = c(16, 17), name = "Model w lowest AIC") +
  geom_errorbar(aes(ymin=M_AIC.2, ymax=M_AIC.3, color=target)) +
  geom_hline(yintercept = 1) +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate") + xlab("Experiment ID")
p_M_antibiotic

# Mutant fitness cost or heterogeneous stress responses explaining the data?
print(df_SIM$by_AIC)

p_S_antibiotic <- ggplot(data = df_SIM, aes(x=ID, y=S_AIC.1, group=target)) + 
  geom_point(aes(color=target, shape = homhet), size=2) + scale_shape_manual(values = c(16, 17), name = "Model w lowest AIC") +
  geom_errorbar(aes(ymin=S_AIC.2, ymax=S_AIC.3, color=target)) +
  geom_hline(yintercept = 1) +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate") + xlab("Experiment ID")
p_S_antibiotic

# Model selection for experiments with significant increase in population-wide mutation rate
selected_models <- data.frame(antibiotic=rep(unique(df_SIM$antibiotic), each=2))
selected_models$target <- mapvalues(selected_models$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target_group))
selected_models$m <- rep(c("hom","het"), length(unique(df_SIM$antibiotic)))
criterion <- "by_AIC_corr"
n <- match(criterion, names(df_SIM))
v <- numeric(length(selected_models$antibiotic))
for (i in 1:length(unique(df_SIM$antibiotic))) {
  v[2*i-1] <- sum(is.element(subset(df_SIM, antibiotic == selected_models$antibiotic[2*i])[,n], c("hom_wo_fitm","hom_fitm","hom_fitm_unconstr"))) 
  v[2*i] <- sum(is.element(subset(df_SIM, antibiotic == selected_models$antibiotic[2*i])[,n], c("het_zero_div","het_div","het_div_fon")))
}
selected_models$prevalence <- v

p_msel_a <- ggplot(data = selected_models, aes(x=factor(m, c("hom","het")), y=prevalence, fill=antibiotic)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial") +
  xlab("Selected model") + ylab("Number of experiments") + scale_x_discrete(labels = c("homogeneous", "heterogeneous"))
p_msel_a

p_msel_t <- ggplot(data = selected_models, aes(x=factor(m, c("hom","het")), y=prevalence, fill=target)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = turbo(length(unique(df_SIM$target)))) +
  xlab("Selected model") + ggtitle(criterion) + ylab("Number of experiments") 
p_msel_t

# Difference in AIC_corr between homogeneous and heterogeneous-response model
p_Delta <- ggplot(data = df_SIM, aes(x=ID, y=Delta_AIC_corr, group=antibiotic)) + geom_point(aes(color=antibiotic, shape=homhet)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  geom_hline(yintercept = 2, linetype = "dashed") + geom_hline(yintercept = -2, linetype = "dashed") + 
  ylab("Difference in AIC") + xlab("Experiment ID") + scale_shape_manual(values = c(16, 17), name = "Model w lowest AIC") 
p_Delta

# Case studies
i <- 20
id <- df_SIM$ID[i]
print(c(as.character(id), df_SIM$by_AIC[i]))
df_id <- read.csv(paste0("experimental_data/model_selection/est_sum_", id, ".csv"))[,-1]
print(as.character(subset(df_SIM, by_AIC == "hom_fitm_unconstr")$ID))
print(as.character(subset(df_SIM, by_AIC == "hom_fitm")$ID))

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
  geom_point(data = df_c, aes(x=concentration, y=M_wo_fitm.1)) +
  geom_errorbar(data = df_c, aes(x=concentration, ymin=M_wo_fitm.2, ymax=M_wo_fitm.3)) +
  geom_hline(yintercept = 1) + scale_y_continuous(trans="log10") +
  geom_point(data = df_Pribis, aes(x=cipro_concentration/1000, y=fold_induction_mutation_rate), color='blue') +
  geom_errorbar(data = df_Pribis, aes(x=cipro_concentration/1000, ymin=fold_min, ymax=fold_max), color='blue') +
  ylab("Fold-change population-wide mutation rate") + ggtitle("Cip")
p_c
