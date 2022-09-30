###################################################
###          Install required packages          ###
###################################################
library(dplyr)
library(ggplot2)
library(ggsignif)

###################################################
###                 Import Data                 ###
###################################################

#setwd()
data_merge_bi <- read.csv("databases/fig4ac.csv")
data_goiK13 <- read.csv("databases/fig4bd.csv")

############
# Figure 4 # run 
############
# CRN WT v Mutant: N, median, P-values
data_merge_bi %>%
  group_by(any_coronin) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE),
            age_q25 = quantile(weighted_RSA,probs = .25,  na.rm = TRUE),
            age_q75 = quantile(weighted_RSA,probs = .75,  na.rm = TRUE))
wilcox.test(weighted_RSA ~ any_coronin, data_merge_bi, conf.int = TRUE)

# CRN WT v Mutant: Figure
crn <-
  ggplot(subset(data_merge_bi, !is.na(any_coronin)), aes(as.factor(any_coronin), weighted_RSA, color = as.factor(any_coronin))) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("CRN Genotype") +
  geom_signif(y_position = 47, xmin = 1, xmax = 2,
              annotation = c("p = 0.04"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","Mut")) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold")) +
  labs(tag = "A")
crn  
ggsave(plot = crn, "figures/Fig4a_crnvRSA.tiff", width = 4, height = 3.5 )
ggsave(plot = crn, "figures/Fig4a_crnvRSA.eps", width = 4, height = 3.5 )

# CRN+K13: N, median, P-values
data_goiK13 %>%
  group_by(any_coronin) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE),
            age_q25 = quantile(weighted_RSA, probs = .25,  na.rm = TRUE),
            age_q75 = quantile(weighted_RSA, probs = .75,  na.rm = TRUE))

wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==0|any_coronin==1))
wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==0|any_coronin==2))
wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==0|any_coronin==3))
wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==1|any_coronin==2))
wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==1|any_coronin==3))
wilcox.test(weighted_RSA ~ any_coronin, subset(data_goiK13, any_coronin==2|any_coronin==3))

# CRN+K13: Figure
crn469 <-
  ggplot(subset(data_goiK13, !is.na(any_coronin) & (any_coronin<5)), aes(as.factor(any_coronin), weighted_RSA, color=as.factor(any_coronin))) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("Genotype") +
  geom_signif(y_position = 47, xmin = 1, xmax = 4,
              annotation = c("p = 0.009"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 57, xmin = 3, xmax = 4,
              annotation = c("p = 0.018"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("C469+CRN-WT","C469+CRN-Mut", "469Y+CRN-WT", "469Y+CRN-Mut")) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold"))+
  labs(tag = "B")
crn469  
ggsave(plot = crn469, "figures/Fig4b_crn469vRSA.tiff", width = 6, height = 3.5 )
ggsave(plot = crn469, "figures/Fig4b_crn469vRSA.eps", width = 6, height = 3.5 )

# CRN WT v Mutant: N, median, P-values
data_merge_bi %>%
  group_by(any_fp2a) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE),
            age_q25 = quantile(weighted_RSA,probs = .25,  na.rm = TRUE),
            age_q75 = quantile(weighted_RSA,probs = .75,  na.rm = TRUE))
wilcox.test(weighted_RSA ~ any_fp2a, data_merge_bi, conf.int = TRUE)

# FP2a WT v Mutant: Figure
fp2a <-
  ggplot(subset(data_merge_bi, !is.na(any_fp2a)), aes(as.factor(any_fp2a), weighted_RSA, color = as.factor(any_fp2a))) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("FP2a Genotype") +
  geom_signif(y_position = 47, xmin = 1, xmax = 2,
              annotation = c("p = 0.015"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("WT","Mut")) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold"))+
  labs(tag = "C")
fp2a  
ggsave(plot = fp2a, "figures/Fig4c_FP2avRSA.tiff", width = 4, height = 3.5 )
ggsave(plot = fp2a, "figures/Fig4c_FP2avRSA.eps", width = 4, height = 3.5 )

# CRN+K13: N, median, P-values
data_goiK13 %>%
  group_by(any_fp2a) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE),
            age_q25 = quantile(weighted_RSA, probs = .25,  na.rm = TRUE),
            age_q75 = quantile(weighted_RSA, probs = .75,  na.rm = TRUE))

wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==0|any_fp2a==1))
wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==0|any_fp2a==2))
wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==0|any_fp2a==3))
wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==1|any_fp2a==2))
wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==1|any_fp2a==3))
wilcox.test(weighted_RSA ~ any_fp2a, subset(data_goiK13, any_fp2a==2|any_fp2a==3))

# CRN+K13: Figure
fp2a469 <-
  ggplot(subset(data_goiK13, !is.na(any_fp2a) & any_fp2a < 4 & any_fp2a != "NA"), aes(as.factor(any_fp2a), weighted_RSA, color=as.factor(any_fp2a))) +
  geom_jitter(width = 0.1) +
  geom_boxplot(outlier.shape = NA, alpha = .2) +
  ylab("RSA Survival (%)") +
  xlab("Genotype") +
  geom_signif(y_position = 47, xmin = 1, xmax = 2,
              annotation = c("p = 0.007"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 52, xmin = 1, xmax = 3,
              annotation = c("p = 0.015"), textsize = 3, tip_length = .05, color = "black") +
  geom_signif(y_position = 57, xmin = 1, xmax = 4,
              annotation = c("p = 0.004"), textsize = 3, tip_length = .05, color = "black") +
  scale_x_discrete(labels=c("C469+FP2a-WT","C469+FP2a-Mut", "469Y+FP2a-WT", "469Y+FP2a-Mut")) +
  scale_color_brewer(palette = "Set1")+
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(color = c("#e41a1c", "#377eb9", "#4daf4a", "#984ea3", "#ff7f00"), face="bold"))+
  labs(tag = "D")
fp2a469  
ggsave(plot = fp2a469, "figures/Fig4d_fp2a469vRSA.tiff", width = 6, height = 3.5 )
ggsave(plot = fp2a469, "figures/Fig4d_fp2a469vRSA.eps", width = 6, height = 3.5 )
