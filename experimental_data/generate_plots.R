library("plyr")
library("ggplot2")
library("ggpubr")

# Read data frames
antibiotic_classes <- read.csv("experimental_data/antibiotic_classes.csv")[,-1]
meta_data <- read.csv("experimental_data/meta_data.csv")[,-1]
sum_data <- read.csv("experimental_data/sum_data.csv")[,-1]
est_paras <- read.csv("experimental_data/est_paras.csv")[,-1]
est_sum <- read.csv("experimental_data/est_sum.csv")[,-1]

# Group mode of action into 'Wall', 'Pro', 'DNA' and 'Other', pool replicates, and order data frames
meta_data$mode_of_action <- mapvalues(meta_data$antibiotic, from = antibiotic_classes$antibiotic_abbr, to = antibiotic_classes$my_group)
meta_data <- subset(meta_data, replicate == 0)
df <- merge(meta_data, est_sum, by = "ID")
df <- arrange(df, mode_of_action, antibiotic)
df$antibiotic <- factor(df$antibiotic, levels = unique(df$antibiotic), ordered = TRUE)
df$ID <- factor(df$ID, levels = unique(df$ID), ordered = TRUE)
df$mode_of_action <- factor(df$mode_of_action, levels = unique(df$mode_of_action), ordered = TRUE)
antibiotic_classes <- arrange(antibiotic_classes, my_group, antibiotic_abbr)
antibiotic_classes$antibiotic_abbr <- factor(antibiotic_classes$antibiotic_abbr, levels = unique(antibiotic_classes$antibiotic_abbr), ordered = TRUE)

# Color-coding antibiotics used in the studies by mode of action 
u <- unique(antibiotic_classes$my_group)
antibiotic_classes$color <- ""
for (g in unique(antibiotic_classes$my_group)){
  a <- antibiotic_classes$antibiotic[antibiotic_classes$my_group == g]
  for (j in 1:length(a)){
    if (g == "DNA"){
      antibiotic_classes$color[antibiotic_classes$antibiotic == a[j]] <- colorRampPalette(c("#00FFFF", "#000112"))(length(a))[j]
    } else {
      if (g == "Wall"){
        antibiotic_classes$color[antibiotic_classes$antibiotic == a[j]] <- colorRampPalette(c("#008631", "#abf7b1"))(length(a))[j]
      } else {
        if (g == "Pro"){
          antibiotic_classes$color[antibiotic_classes$antibiotic == a[j]] <- colorRampPalette(c("#FFFF9E", "#b30000"))(length(a))[j]
        } else {
          antibiotic_classes$color[antibiotic_classes$antibiotic == a[j]] <- colorRampPalette(c("#FF00CC", "#CC00CC"))(length(a))[j]
        }
      }
    }
  }
}
v <- numeric(length(antibiotic_classes$antibiotic))
for (i in 1:length(v)) {
  v[i] <- length(subset(meta_data, antibiotic == antibiotic_classes$antibiotic_abbr[i])$antibiotic)
}
antibiotic_classes$prevalence <- v
p_antibiotic <- ggplot(data = antibiotic_classes, aes(x=my_group, y=prevalence, fill=antibiotic_abbr)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = antibiotic_classes$color) + xlab("Mode of action")
p_antibiotic

# Estimated increase in population-wide mutation rate by antibiotic
df_SIM <- subset(df, is.element(SIM, c("yes", "hom")))
p_M_antibiotic <- ggplot(data = df_SIM, aes(x=ID, y=M.1, group=antibiotic)) + 
  geom_point(aes(color=antibiotic)) + geom_point(aes(color=antibiotic)) +
  geom_errorbar(aes(ymin=M.2, ymax=M.3, color=antibiotic)) +
  geom_hline(yintercept = 1) +
  scale_color_manual(values = subset(antibiotic_classes, is.element(antibiotic_abbr, unique(df_SIM$antibiotic)))$color) + 
  scale_y_continuous(trans="log10") + theme(axis.text.x = element_text(angle = 90), plot.margin = margin(3.5,0.5,0.5,0.5, "cm")) +
  ylab("Fold change population-wide mutation rate")
p_M_antibiotic

# Estimated increase in population-wide mutation rate by mode of action
# Model used in the inference, for purpose of comparison
m <- "hom_wo_fitm"
df_W <- merge(subset(est_paras, model == m), meta_data, by = "ID")

# Testing for normality -> not normal
shapiro.test(subset(df_W, mode_of_action == "DNA")$M_MLE)

# Kruskal-Wallis test and pairwise comparison -> DNA-Pro and DNA-Wall are significantly different
kruskal.test(df_W$M_MLE, df_W$mode_of_action)
pairwise.wilcox.test(df_W$M_MLE, df_W$mode_of_action)
compare_means(M_MLE ~ mode_of_action, data = df_W)
my_comparison <- list(c("DNA","Other"),c("DNA","Pro"),c("DNA","Wall"))
p_M_mode <- ggplot(data = df_W, aes(x=mode_of_action, y=M_MLE, fill=mode_of_action)) + geom_boxplot() +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  scale_fill_manual(values = c("#009092", "#FF00CC", "#FF6600", "#5ced73")) + 
  ylab("Estimated fold change in population-wide mutation rate") + 
  stat_compare_means(comparisons = my_comparison, label.x = 1, label.y = c(70,110,170)) + stat_compare_means(label.y = 400) 
p_M_mode
# Comparing DNA - not DNA
levels(df_W$mode_of_action) <- c("DNA", "Other", "Pro", "Wall", "not DNA")
df_W$mode_of_action[which(is.element(df_W$mode_of_action, c("Pro","Other","Wall")))] <- "not DNA"
wilcox.test(M_MLE ~ mode_of_action, data = df_W)
p_M_DNA <- ggplot(data = df_W, aes(x=mode_of_action, y=M_MLE, fill=mode_of_action)) + geom_boxplot() +
  coord_trans(y = "log10", ylim = c(5*10^-2,5*10^2)) + scale_y_continuous(breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
  scale_fill_manual(values = c("#009092", "#FF00CC", "#FF6600", "#5ced73")) + 
  ylab("Estimated fold change in population-wide mutation rate") +
  stat_compare_means(label.y = 400) 
p_M_DNA

# Model selection for experiments with significant increase in population-wide mutation rate
selected_models <- data.frame(antibiotic=rep(unique(df_SIM$antibiotic), each=3))
selected_models$mode_of_action <- mapvalues(selected_models$antibiotic, from=antibiotic_classes$antibiotic_abbr, to=antibiotic_classes$my_group) 
selected_models$m <- rep(c("hom","none","het"), length(unique(df_SIM$antibiotic)))
criterion <- "by_AIC"
n <- match(criterion, names(df_SIM))
m <- match("SIM", names(df_SIM))
v <- numeric(length(selected_models$antibiotic))
for (i in 1:length(unique(df_SIM$antibiotic))) {
  v[3*i-2] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "hom") + sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,m] == "hom")
  v[3*i-1] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "none") 
  v[3*i] <- sum(subset(df_SIM, antibiotic == selected_models$antibiotic[3*i])[,n] == "het")
}
selected_models$prevalence <- v

p_model_selection <- ggplot(data = selected_models, aes(x=factor(m, c("hom","none","het")), y=prevalence, fill=mode_of_action)) + geom_bar(stat = "identity")+ 
  scale_fill_manual(values = c("#009092", "#FF00CC", "#FF6600", "#5ced73")) + xlab("Selected model") + ggtitle(criterion)
p_model_selection

# Experiments, for which a heterogeneous stress response is selected or homogeneous/heterogeneous response cannot be distinguished clearly
df_het <- subset(df_SIM, is.element(by_AIC, c("none", "het")))

# Frenoy et al. 2018 Norfloxacin
df_Nor <- subset(est_paras, ID == "Frenoy_Nor")
# Untreated condition
mc_data <- read.table(paste0("experimental_data/raw_counts/Frenoy_LB.txt"), header = FALSE, sep = ",", fill = TRUE)
mc_UT <- read_counts(mc_data[2,])
Nf_UT <- mean(read_counts(mc_data[3,]))
eff_UT <- as.numeric(mc_data[4,1])
p_mc_hom_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[4], eff=eff_UT, fit_m=df_Nor$fitm_UT_MLE[4]) * length(mc_UT)
p_mc_hom_constr_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_UT_MLE[3], eff=eff_UT, fit_m=df_Nor$fitm_UT_MLE[3]) * length(mc_UT)
p_mc_het_UT <- pMudi(max(mc_UT), Nf_UT, df_Nor$mu_off_MLE[6], eff=eff_UT) * length(mc_UT)
print(c(df_Nor$mu_UT_MLE[4], df_Nor$fitm_UT_MLE[4], df_Nor$mu_UT_MLE[3], df_Nor$fitm_UT_MLE[3]))
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
p_mc_hom_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[4], eff=eff_S, fit_m=df_Nor$fitm_S_MLE[4]) * length(mc_S)
p_mc_hom_constr_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_S_MLE[3], eff=eff_S, fit_m=df_Nor$fitm_S_MLE[3]) * length(mc_S)
p_mc_het_S <- pMudi(max(mc_S), Nf_S, df_Nor$mu_off_MLE[6], eff=eff_S, S=df_Nor$S_MLE[6], f_on=df_Nor$f_on_MLE[6], rel_div_on=df_Nor$rel_div_on_MLE[6]) * length(mc_S)
print(c(df_Nor$mu_S_MLE[4], df_Nor$fitm_S_MLE[4], df_Nor$mu_S_MLE[3], df_Nor$fitm_S_MLE[3]))
print(c(df_Nor$mu_off_MLE[6], df_Nor$S_MLE[6], df_Nor$f_on_MLE[6], df_Nor$rel_div_on_MLE[6]))
p_mc_s <- ggplot() + geom_bar(aes(mc_S), fill="#89CFF0") + xlab("Number of colonies") + ylab("Number of plates") +
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_S), color="#820298") + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_hom_constr_S), color="#f1aafd") + 
  geom_line(aes(x=0:max(mc_S), y=p_mc_het_S), color="#FFD300") + ggtitle("Norfloxacin")
p_mc_s

# Dose-dependence
# Model used in the inference, for purpose of comparison
m <- "hom_wo_fitm"
df_W <- merge(subset(est_paras, model == m), meta_data, by = "ID")
df_Pribis <- read.delim("experimental_data/Pribis_Fig2D.txt", header = TRUE, sep = '\t', comment.char="#")
df_c <- subset(subset(subset(df_W, microbe == "E. coli"), strain != "DeltaLexA"), antibiotic == "Cip")
p_c <- ggplot() + 
  geom_point(data = df_c, aes(x=of_MIC, y=M_MLE)) +
  geom_errorbar(data = df_c, aes(x=of_MIC, ymin=M_lower_bound, ymax=M_upper_bound)) +
  geom_hline(yintercept = 1) +
  geom_point(data = df_Pribis, aes(x=of_MIC, y=fold_induction_mutation_rate), color='blue') +
  geom_errorbar(data = df_Pribis, aes(x=of_MIC, ymin=fold_min, ymax=fold_max), color='blue') +
  ylab("Fold change population-wide mutation rate") + ggtitle(a)
p_c
