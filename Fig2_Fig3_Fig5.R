###################################################
###          Install required packages          ###
###################################################
library(ggplot2)
library(ggsignif)
library(patchwork)

###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())

#setwd("2022 Tumwebaze et al Nature comm")

data_wide <- read.csv("databases/exvivo_target-k13_MayAug21_wide.csv")


###################################################
###                    Analyses                 ###
###################################################
## Figure 2:  LM and RSA by region
##########################
# N, median, p-values
subset(data_wide, !is.na(LM)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))
wilcox.test(LM ~ site3, data_wide, exact = FALSE, conf.int = TRUE)

LMvregion <-
  ggplot(subset(data_wide, !is.na(LM)), aes(site3, LM, color=as.factor(site3))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("Lumefantrine IC50 (nM)") +
  xlab("Region") +
  geom_signif(y_position = 97, xmin = 1, xmax = 2,
              annotation = c("p < 0.0001"), textsize = 3, tip_length = .05, color = "black") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9"), face="bold")) +
  labs(tag = "A")
LMvregion
ggsave(plot = LMvregion, "figures/Fig2_LMvRegion.tiff", width = 3.5, height = 3.5 )

# N, median, p-values
subset(data_wide, !is.na(weighted_RSA) & include==1) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ site3, subset(data_wide, include ==1), conf.int = TRUE)

RSAvregion <-
  ggplot(subset(data_wide, !is.na(weighted_RSA) & include==1), aes(site3, weighted_RSA, color=as.factor(site3))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("Region") +
  geom_signif(y_position = 50, xmin = 1, xmax = 2,
              annotation = c("p = 0.53"), textsize = 3, tip_length = .05, color = "black") +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9"), face="bold")) +
  labs(tag = "B")
RSAvregion
ggsave(plot = RSAvregion, "figures/Fig2_RSAvRegion.tiff", width = 3.5, height = 3.5 )


##########################
## Figure 3: RSA by genotype
##########################
# N, median, P-values; combined mix and mutant
subset(data_wide, include ==1) %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ k13_genotype, subset(data_wide, include == 1 & (k13_genotype == 0 | k13_genotype==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype, subset(data_wide, include == 1 & (k13_genotype == 0 | k13_genotype==2)), conf.int = TRUE)

RSAvK13genotype <-
  ggplot(subset(data_wide, !is.na(k13_genotype) & !is.na(weighted_RSA) & include ==1), aes(as.factor(k13_genotype), weighted_RSA, color = as.factor(k13_genotype))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("K13 Genotype") +
  geom_signif(y_position = 50, xmin = 1, xmax = 2,
              annotation = c("p = 0.03"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 55, xmin = 1, xmax = 3,
              annotation = c("p = 0.35"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","C469Y", "A675V")) +
  scale_color_manual(values = c("#e41a1c", "#377eb9", "#984ea3")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#984ea3"), face="bold"))+
  labs(tag = "A")
RSAvK13genotype  
ggsave(plot = RSAvK13genotype, "figures/Fig3a_RSAvK13genotype.tiff", width = 4, height = 3.5 )

# N, median, P-values; mix and mutant separate
subset(data_wide, include ==1) %>%
  group_by(k13_genotype2) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==2)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==3)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==4)), conf.int = TRUE)

RSAvK13genotype2 <-
  ggplot(subset(data_wide, !is.na(k13_genotype2) & !is.na(weighted_RSA) & include==1), aes(as.factor(k13_genotype2), weighted_RSA, color = as.factor(k13_genotype2))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("K13 Genotype") +
  geom_signif(y_position = 50, xmin = 1, xmax = 2,
              annotation = c("p = 0.47"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 55, xmin = 1, xmax = 3,
              annotation = c("p = 0.009"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 60, xmin = 1, xmax = 4,
              annotation = c("p = 0.37"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 65, xmin = 1, xmax = 5,
              annotation = c("p = 0.82"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","C469Y", "469Y","A675V", "675V")) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold")) +
  labs(tag = "B")
RSAvK13genotype2
ggsave(plot = RSAvK13genotype2, "figures/Fig3b_RSAvK13genotype.tiff", width = 4, height = 3.5 )




##########################
## Figure 5: Lumefantrine by genotype
##########################
# N, median, P-values; combined mix and mutant
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))
wilcox.test(LM ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(LM ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

LMvK13genotype <-
  ggplot(subset(data_wide, !is.na(k13_genotype) & !is.na(LM)), aes(as.factor(k13_genotype), LM, color=as.factor(k13_genotype))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("Lumefantrine IC50 (nM)") +
  xlab("K13 Genotype") +
  geom_signif(y_position = 83, xmin = 1, xmax = 2,
              annotation = c("p = 0.02"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 90, xmin = 1, xmax = 3,
              annotation = c("p = 0.02"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","C469Y", "A675V")) +
  scale_color_manual(values = c("#e41a1c", "#377eb9", "#984ea3")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#984ea3"), face="bold"))+
  labs(tag = "A")

LMvK13genotype  
ggsave(plot = LMvK13genotype, "figures/Fig5a_LMvK13genotype.tiff", width = 4, height = 3.5 )

# N, median, P-values
subset(data_wide) %>%
  group_by(k13_genotype2) %>%
  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))
wilcox.test(LM ~ k13_genotype2, subset(data_wide, (k13_genotype2 == 0 | k13_genotype2==1)), conf.int = TRUE)
wilcox.test(LM ~ k13_genotype2, subset(data_wide, (k13_genotype2 == 0 | k13_genotype2==2)), conf.int = TRUE)
wilcox.test(LM ~ k13_genotype2, subset(data_wide, (k13_genotype2 == 0 | k13_genotype2==3)), conf.int = TRUE)
wilcox.test(LM ~ k13_genotype2, subset(data_wide, (k13_genotype2 == 0 | k13_genotype2==4)), conf.int = TRUE)

LMvK13genotype2 <-
  ggplot(subset(data_wide, !is.na(k13_genotype2) & !is.na(LM)), aes(as.factor(k13_genotype2), LM, color=as.factor(k13_genotype2))) +
  geom_jitter() +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("Lumefantrine IC50 (nM)") +
  xlab("K13 Genotype") +
  geom_signif(y_position = 82, xmin = 1, xmax = 2,
              annotation = c("p = 0.21"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 89, xmin = 1, xmax = 3,
              annotation = c("p = 0.02"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 96, xmin = 1, xmax = 4,
              annotation = c("p = 0.04"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 103, xmin = 1, xmax = 5,
              annotation = c("p = 0.27"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","C469Y", "469Y","A675V", "675V"), ) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold")) +
  labs(tag = "B")
LMvK13genotype2  
ggsave(plot = LMvK13genotype2, "figures/Fig5b_LMvK13genotype.tiff", width = 4, height = 3.5 )




