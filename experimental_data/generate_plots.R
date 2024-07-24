library("plyr")
library("ggplot2")
library("viridisLite")
library("ggpubr")

# Read data frames
antibiotic_classes <- read.csv("experimental_data/antibiotic_classes.csv")[,-1]
meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]
sum_data <- read.csv("experimental_data/sum_data.csv")[,-1]
est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
est_sum <- read.csv("experimental_data/est_sum.csv")[,-1]

# Add target to meta data, pool replicates, and order data frames
meta_data$target <- mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$target_group)
meta_data <- subset(meta_data, replicate == 0)
df <- merge(meta_data, sum_data, by = "ID")
df <- merge(df, est_sum, by = "ID")
df$antibiotic <- factor(df$antibiotic, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)
df$target <- factor(df$target, levels = unique(antibiotic_classes$target_group), ordered = TRUE)
antibiotic_classes$antibiotic_abbr <- factor(antibiotic_classes$antibiotic_abbr, levels = antibiotic_classes$antibiotic_abbr, ordered = TRUE)
antibiotic_classes$target_group <- factor(antibiotic_classes$target_group, levels = unique(antibiotic_classes$target_group), ordered = TRUE)

# Color-coding antibiotics used in the studies (by target group)
antibiotic_classes$color <- turbo(length(antibiotic_classes$antibiotic_abbr))
v <- numeric(length(antibiotic_classes$antibiotic_abbr))
for (i in 1:length(v)) {
  v[i] <- length(subset(meta_data, antibiotic == antibiotic_classes$antibiotic_abbr[i])$antibiotic)
}
antibiotic_classes$prevalence <- v
p_antibiotic <- ggplot(data = antibiotic_classes, aes(x=target_group, y=prevalence, fill=antibiotic_abbr)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = antibiotic_classes$color) + xlab("Grouped by target") + ylab("Number of experiments")
p_antibiotic

# Estimated increase in population-wide mutation rate by antibiotic
# Homogeneous-response model without differential mutant fitness used in the inference, for purpose of comparison
p_M_antibiotic <- ggplot(data = df, aes(x=ID, y=M_wo_fitm.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M_wo_fitm.2, ymax=M_wo_fitm.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df$antibiotic)))$color) + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold-change population-wide mutation rate")
p_M_antibiotic

# Experiments for which SIM was detected
df_SIM <- subset(df, SIM == TRUE)
p_M_antibiotic <- ggplot(data = df_SIM, aes(x=ID, y=M.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M.2, ymax=M.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color) + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Increase population-wide mutation rate")
p_M_antibiotic

# Testing for normality -> not normal
shapiro.test(subset(df, is.element(, c("DNA", "Gyrase")))$M_wo_fitm.1)

# Kruskal-Wallis test and pairwise comparison -> DNA and ribosome binding significantly different
kruskal.test(df$M_wo_fitm.1, df$target)
print(compare_means(M_wo_fitm.1 ~ target, data = df), n = 45)
# Comparing antimicrobials that directly target DNA/gyrase with others
df$DNA_direct_target <- logical(length(df$target))
df$DNA_direct_target[is.element(df$target, c("DNA", "Gyrase"))] <- TRUE
kruskal.test(df$M_wo_fitm.1, df$DNA_direct_target)
wilcox.test(M_wo_fitm.1 ~ DNA_direct_target, data = df)
p_M_DNA <- ggplot(data = df, aes(x=DNA_direct_target, y=M_wo_fitm.1, fill=DNA_direct_target)) + geom_boxplot() +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  ylab("Estimated fold change in population-wide mutation rate") +
  stat_compare_means(label.y = 400) 
p_M_DNA

# Model selection for experiments with significant increase in population-wide mutation rate
selected_models <- data.frame(antibiotic=rep(unique(df_SIM$antibiotic), each=3))
selected_models$target <- mapvalues(selected_models$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = as.character(antibiotic_classes$target_group))
selected_models$m <- rep(c("hom","none","het"), length(unique(df_SIM$antibiotic)))
criterion <- "by_AIC"
n <- match(criterion, names(df_SIM))
v <- numeric(length(selected_models$antibiotic))
for (i in 1:length(unique(df_SIM$antibiotic))) {
  v[3*i-2] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "hom") 
  v[3*i-1] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "none") 
  v[3*i] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "het")
}
selected_models$prevalence <- v

p_msel_a <- ggplot(data = selected_models, aes(x=factor(m, c("hom","none","het")), y=prevalence, fill=antibiotic)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color) +
  xlab("Selected model") + ggtitle(criterion) + ylab("Number of experiments")
p_msel_a

p_msel_t <- ggplot(data = selected_models, aes(x=factor(m, c("hom","none","het")), y=prevalence, fill=target)) + geom_bar(stat = "identity") + 
  scale_fill_manual(values = turbo(length(unique(df_SIM$target)))) +
  xlab("Selected model") + ggtitle(criterion) + ylab("Number of experiments") 
p_msel_t

# Difference in AIC between homogeneous and heterogeneous-response model
p_Delta_AIC <- ggplot(data = df_SIM, aes(x=ID, y=Delta_AIC, group=antibiotic)) + geom_point(aes(color=antibiotic)) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color) +
  theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  geom_hline(yintercept = 2) + geom_hline(yintercept = -2)
p_Delta_AIC

# Experiments, for which a heterogeneous stress response is selected or homogeneous/heterogeneous response cannot be distinguished clearly
df_het <- subset(df_SIM, is.element(by_AIC, c("none", "het")))

# Frenoy et al. 2018 Norfloxacin
df_Nor <- subset(est_paras, ID == "Frenoy_Nor")
print(subset(est_sum, ID == "Frenoy_Nor"))
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_LB.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])
Nf_UT <- mean(read_counts(mc_data[3,]))
eff_UT <- as.numeric(mc_data[4,1])
p_mc_hom_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[5], eff=eff_UT, fit_m=df_Nor$fitm_UT_MLE[5]) * length(mc_UT)
p_mc_hom_constr_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[3], eff=eff_UT) * length(mc_UT)
p_mc_het_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_off_MLE[7], eff=eff_UT) * length(mc_UT)
print(c(df_Nor$mu_UT_MLE[5], df_Nor$fitm_UT_MLE[5], df_Nor$mu_UT_MLE[3]))
print(df_Nor$mu_off_MLE[6])
p_mc_UT <- ggplot() + geom_bar(aes(mc_UT), fill="#123288") + xlab("Number of colonies") + ylab("Number of plates") +
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_UT), color="#820298") + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_hom_constr_UT), color="#f1aafd") + 
  geom_line(aes(x=0:max(mc_UT), y=p_mc_het_UT), color="#FFD300") + ggtitle("Untreated condition")
p_mc_UT
# Stressful condition: 0.05 mug/mL Norfloxacin
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_Nor.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_S <- read_counts(mc_data[2,])
Nf_S <- mean(read_counts(mc_data[3,]))
eff_S <- as.numeric(mc_data[4,1])
p_mc_hom_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[5], eff=eff_S, fit_m=df_Nor$fitm_S_MLE[5]) * length(mc_S)
p_mc_hom_constr_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[3], eff=eff_S) * length(mc_S)
p_mc_het_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_off_MLE[7], eff=eff_S, S=df_Nor$S_MLE[7], f_on=df_Nor$f_on_MLE[7], rel_div_on=df_Nor$rel_div_on_MLE[7]) * length(mc_S)
print(c(df_Nor$mu_S_MLE[5], df_Nor$fitm_S_MLE[5], df_Nor$mu_S_MLE[3]))
print(c(df_Nor$mu_off_MLE[7], df_Nor$S_MLE[7], df_Nor$f_on_MLE[7], df_Nor$rel_div_on_MLE[7]))
p_mc_s <- ggplot() + geom_bar(aes(mc_S), fill="#89CFF0") + xlab("Number of colonies") + ylab("Number of plates") +
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_S), color="#820298") + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_constr_S), color="#f1aafd") + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_het_S), color="#FFD300") + ggtitle("Norfloxacin")
p_mc_s

# Dose-dependence
df_Pribis <- read.delim("experimental_data/Pribis_Fig2D.txt", header = TRUE, sep = '\t', comment.char="#")
df_c <- subset(subset(subset(df, microbe == "E. coli"), strain != "DeltaLexA"), antibiotic == "Cip")
p_c <- ggplot() + 
  geom_point(data = df_c, aes(x=concentration, y=M_wo_fitm.1)) +
  geom_errorbar(data = df_c, aes(x=concentration, ymin=M_wo_fitm.2, ymax=M_wo_fitm.3)) +
  geom_hline(yintercept = 1) + scale_y_continuous(trans="log10") +
  geom_point(data = df_Pribis, aes(x=cipro_concentration/1000, y=fold_induction_mutation_rate), color='blue') +
  geom_errorbar(data = df_Pribis, aes(x=cipro_concentration/1000, ymin=fold_min, ymax=fold_max), color='blue') +
  ylab("Fold-change population-wide mutation rate") + ggtitle("Cip")
p_c
