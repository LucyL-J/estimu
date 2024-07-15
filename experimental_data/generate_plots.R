library("plyr")
library("ggplot2")

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
