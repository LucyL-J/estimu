library("plyr")
library("ggplot2")
library("viridisLite")
library("ggpubr")
library("lme4")
library("performance")

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
df <- arrange(df, species)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)

# Estimated increase in population-wide mutation rate by antibiotic
# Homogeneous-response model without differential mutant fitness used in the inference, for purpose of comparison
p_M_antibiotic <- ggplot(data = df, aes(x=ID, y=M_wo_fitm.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_wo_fitm.2, ymax=M_wo_fitm.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

df_sel <- df
df_sel$M.1[is.na(df_sel$M.1)] <- 1.
df_sel$M.2[is.na(df_sel$M.2)] <- 1.
df_sel$M.3[is.na(df_sel$M.3)] <- 1.
p_M_antibiotic <- ggplot(data = df_sel, aes(x=ID, y=M.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M.2, ymax=M.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color, name = "Antibiotic (abbr)") + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate") + xlab("Experiment ID") +
  theme(axis.text.x = element_text(vjust = 0.5))
p_M_antibiotic

# Plating efficiency/number of parallel cultures and width of confidence intervals

length(subset(df, plated_fraction < 1)$ID)
df$width_CI <- (df$M_wo_fitm.3-df$M_wo_fitm.2)/df$M_wo_fitm.1
reg <- lm(formula = log10(width_CI) ~ log10(plated_fraction*n_cultures_tot), data = df)
summary(reg)
mlm <- lmer(log10(width_CI) ~ log10(plated_fraction) + log10(n_cultures_tot) + (1|author) + (1|strain), data = df)
summary(mlm)
check_model(mlm)
p_CI <- ggplot(data = df, aes(x=log10(plated_fraction*n_cultures_tot), y=log10(width_CI))) + 
  geom_point() + geom_smooth(method = lm)
p_CI

cor(df[, c("width_CI", "plated_fraction", "n_cultures_tot")], method = "kendall")
cor.test(df$width_CI, df$plated_fraction, method = "kendall")
cor.test(df$width_CI, df$n_cultures_tot, method = "kendall")
cor.test(df$plated_fraction, df$n_cultures_tot, method = "kendall")
p_CI_corr_p <- ggplot(data = df, aes(x=plated_fraction, y=width_CI)) + 
  geom_point(aes(color=log10(n_cultures_tot))) + scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans="log10") + stat_cor(label.y.npc = "bottom", method = "kendall") +
  labs(x="Plated fraction", y="Normalised width of 95% CI around MLE estimate", color="log(c)")
p_CI_corr_p

p_CI_corr_n <- ggplot(data = df, aes(x=n_cultures_tot, y=width_CI)) + 
  geom_point(aes(color=log10(plated_fraction))) + 
  scale_y_continuous(trans="log10") + scale_x_continuous(trans = "log10") + stat_cor(label.y.npc = "top", label.x.npc = "centre", method = "kendall") +
  labs(x="Total number of parallel cultures", y="Normalised width of 95% CI around MLE estimate", color="log(E)")
p_CI_corr_n

p_CI_M <- ggplot(data = df, aes(x=M_wo_fitm.1, y=width_CI)) + 
  geom_point(aes(color=log10(plated_fraction*n_cultures_tot))) + 
  scale_y_continuous(trans="log10") + scale_x_continuous(trans = "log10") +
  labs(x="Fold-change population-wide mutation rate", y="Normalised width of 95% CI around MLE estimate", color="log(E*c)")
p_CI_M

df_glmm <- subset(df, !is.na(of_MIC))
glmm <- glmer(SIM ~ of_MIC + plated_fraction + n_cultures_tot + target + (1|species), data = df, family=binomial)
check_singularity(glmm, tolerance=10^-7)
summary(glmm)
check_model(glmm)
model_performance(glmm)

my_glm <- glm(SIM ~ of_MIC + plated_fraction + n_cultures_tot + target, data = df, family=binomial)
summary(my_glm)
check_model(my_glm)

print(subset(df, SIM == TRUE)$ID)
print(c(length(subset(df, SIM == TRUE)$ID), print(length(subset(df, M_wo_fitm.1 > 1)$ID))))

# Further analysis with experiments using E. coli MG1655 and TD2158 (no mutant strains)
df <- subset(subset(df, species == "E. coli"), is.element(strain, c("MG1655", "TD2158")))
print(subset(df, SIM == TRUE)$ID)

# Testing for normality -> not normal
shapiro.test(df$M_wo_fitm.1)
length(subset(df, is.element(target, c("DNA", "DNA gyrase")))$ID)
shapiro.test(subset(df, is.element(target, c("DNA", "DNA gyrase")))$M_wo_fitm.1)
length(subset(df, target=="Ribosome")$ID)
shapiro.test(subset(df, target == "Ribosome")$M_wo_fitm.1)

# Kruskal-Wallis test -> DNA/DNA-gyrase and ribosome binding significantly different
kruskal.test(M_wo_fitm.1 ~ target, data = df)
df_KW <- subset(df, is.element(target, c("DNA", "DNA gyrase", "Ribosome")))
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
df_SIM <- subset(df, SIM == TRUE)
df_SIM <- arrange(df_SIM, target)
df_SIM$ID <- factor(df_SIM$ID, levels = unique(df_SIM$ID), ordered = TRUE)
print(c(length(df_SIM$ID), length(subset(df, M_wo_fitm.1>1)$ID)))
print(c(length(subset(df_SIM, is.element(target, c("DNA", "DNA gyrase")))$ID),length(subset(df_SIM, target=="Ribosome")$ID)))
p_M_antibiotic <- ggplot(data = df_SIM, aes(x=ID, y=M.1, group=target)) + 
  geom_point(aes(color=target)) +
  geom_errorbar(aes(ymin=M.2, ymax=M.3, color=target)) +
  geom_hline(yintercept = 1) +
  scale_colour_manual(values = turbo(length(unique(df_SIM$target))), name = "Target") +
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate") + xlab("Experiment ID")
p_M_antibiotic

# Exclude experiments with less than 3 parallel cultures in the fluctuation assay under antimicrobial treatment
df_SIM <- subset(df_SIM, n_cultures >= 3)
length(df_SIM$ID)
print(c(length(subset(df, M_wo_fitm.1>1)$ID), length(df_SIM$ID)))
print(c(length(subset(subset(df, M_wo_fitm.1>1), is.element(target, c("DNA", "Gyrase")))$ID),length(subset(subset(df, M_wo_fitm.1>1), target=="Ribosome")$ID)))
print(c(length(subset(df_SIM, is.element(target, c("DNA", "Gyrase")))$ID),length(subset(df_SIM, target=="Ribosome")$ID)))

# Model selection for experiments with significant increase in population-wide mutation rate
selected_models <- data.frame(antibiotic=rep(unique(df_SIM$antibiotic), each=3))
selected_models$target <- mapvalues(selected_models$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target_group))
selected_models$m <- rep(c("hom","none","het"), length(unique(df_SIM$antibiotic)))
criterion <- "by_BIC"
n <- match(criterion, names(df_SIM))
v <- numeric(length(selected_models$antibiotic))
for (i in 1:length(unique(df_SIM$antibiotic))) {
  v[3*i-2] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "hom") 
  v[3*i-1] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "none") 
  v[3*i] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "het")
}
selected_models$prevalence <- v

p_msel_a <- ggplot(data = selected_models, aes(x=factor(m, c("hom","none","het")), y=prevalence, fill=antibiotic)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial") +
  xlab("Selected model") + ylab("Number of experiments") + scale_x_discrete(labels = c("homogeneous", "none", "heterogeneous"))
p_msel_a

p_msel_t <- ggplot(data = selected_models, aes(x=factor(m, c("hom","none","het")), y=prevalence, fill=target)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = turbo(length(unique(df_SIM$target)))) +
  xlab("Selected model") + ggtitle(criterion) + ylab("Number of experiments") 
p_msel_t

# Difference in AIC between homogeneous and heterogeneous-response model
p_Delta_AIC <- ggplot(data = df_SIM, aes(x=ID, y=Delta_BIC, group=antibiotic)) + geom_point(aes(color=antibiotic)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color, name = "Antimicrobial") +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.9, hjust = 0.9), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  geom_hline(yintercept = 6, linetype = "dashed") + geom_hline(yintercept = -6, linetype = "dashed") +
  ylab("Difference in BIC") + xlab("Experiment ID")
p_Delta_AIC

# Experiments, for which a heterogeneous stress response is selected or homogeneous/heterogeneous response cannot be distinguished clearly
df_het <- subset(df_SIM, is.element(by_AIC, c("het")))
df_none <- subset(df_SIM, is.element(by_AIC, c("none")))
df_hom <- subset(df_SIM, is.element(by_AIC, c("hom")))

# Frenoy et al. 2018 Norfloxacin
df_Nor <- subset(est_paras, ID == "Frenoy_Nor")
print(subset(est_sum, ID == "Frenoy_Nor"))
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_LB.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])
Nf_UT <- mean(read_counts(mc_data[3,]))
eff_UT <- as.numeric(mc_data[4,1])
p_mc_hom_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[5], plateff=eff_UT, fit_m=df_Nor$fitm_UT_MLE[5]) * length(mc_UT)
p_mc_hom_constr_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[3], plateff=eff_UT) * length(mc_UT)
p_mc_het_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_off_MLE[7], plateff=eff_UT) * length(mc_UT)
print(c(df_Nor$mu_UT_MLE[5], df_Nor$fitm_UT_MLE[5], df_Nor$mu_UT_MLE[3]))
print(df_Nor$mu_off_MLE[7])
p_mc_UT <- ggplot() + geom_histogram(aes(mc_UT, y=..density.. * length(mc_UT)), fill="#1C458A", bins = 50) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_UT), color="#8F3F8C", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_constr_UT), color="#E0B0FF", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_het_UT), color="#FFD300", linewidth = 1., linetype = "dashed") + 
  ggtitle("Untreated experiments") + xlab("Number of colonies") + ylab("Number of plates")
p_mc_UT
# Stressful condition: 0.05 mug/mL Norfloxacin
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_Nor.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])
Nf_S <- mean(read_counts(mc_data[3,]))
eff_S <- as.numeric(mc_data[4,1])
p_mc_hom_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[5], plateff=eff_S, fit_m=df_Nor$fitm_S_MLE[5]) * length(mc_S)
p_mc_hom_constr_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[3], plateff=eff_S) * length(mc_S)
p_mc_het_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_off_MLE[7], plateff=eff_S, S=df_Nor$S_MLE[7], f_on=df_Nor$f_on_MLE[7], rel_div_on=df_Nor$rel_div_on_MLE[7]) * length(mc_S)
print(c(df_Nor$mu_S_MLE[5], df_Nor$fitm_S_MLE[5], df_Nor$mu_S_MLE[3]))
print(c(df_Nor$mu_off_MLE[7], df_Nor$S_MLE[7], df_Nor$f_on_MLE[7], df_Nor$rel_div_on_MLE[7]))
p_mc_s <- ggplot() + geom_histogram(aes(mc_S, y=..density.. *length(mc_S)), fill="#8EC44F", bins = 20) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_S), color="#8F3F8C", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_constr_S), color="#E0B0FF", linewidth = 1.) + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_het_S), color="#FFD300", linewidth = 1., linetype = "dashed") +
  xlab("Number of colonies") + ylab("Number of plates") + ggtitle("Experiments with norfloxacin")
p_mc_s

# Dose-dependence
df_Pribis <- read.delim("experimental_data/Pribis_Fig2D.txt", header = TRUE, sep = '\t', comment.char="#")
df_c <- subset(df, antibiotic == "Cip")
p_c <- ggplot() + 
  geom_point(data = df_c, aes(x=concentration, y=M_wo_fitm.1)) +
  geom_errorbar(data = df_c, aes(x=concentration, ymin=M_wo_fitm.2, ymax=M_wo_fitm.3)) +
  geom_hline(yintercept = 1) + scale_y_continuous(trans="log10") +
  geom_point(data = df_Pribis, aes(x=cipro_concentration/1000, y=fold_induction_mutation_rate), color='blue') +
  geom_errorbar(data = df_Pribis, aes(x=cipro_concentration/1000, ymin=fold_min, ymax=fold_max), color='blue') +
  ylab("Fold-change population-wide mutation rate") + ggtitle("Cip")
p_c
