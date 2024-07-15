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