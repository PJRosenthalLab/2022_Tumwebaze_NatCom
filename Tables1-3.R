

###################################################
###          Install required packages          ###
###################################################
library(dplyr)
library(psych)
library(stringr)
library(naniar)
library(tidyr)
library(boot)
library(lubridate)
library(ggplot2)
###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())

#setwd("")

data_wide <- read.csv("databases/exvivo_target-k13_MayAug21_wide.csv")

data_wide$k13_genotype <- ifelse(data_wide$all_469==1 | data_wide$all_469==2, 1, 
                                 ifelse(data_wide$all_675V ==1 | data_wide$all_675V ==2, 2, 
                                        ifelse(data_wide$all_WT_469==0, 0, NA)))
data_wide$k13_genotype <- ifelse(data_wide$all_469==1 & data_wide$all_675V==1, 2, data_wide$k13_genotype)

data_wide$all_WT3_469 <- ifelse(data_wide$all_469Y==2, 2, data_wide$all_WT_469)
data_wide$all_WT3_675 <- ifelse(data_wide$all_675V==2, 2, data_wide$all_WT_675)

subset(data_wide, (all_469Y==1|all_469Y==2) & (all_675V==1|all_675V==2), select = c("id", "all_469Y", "all_675V"))

data_wide$k13_genotype2 <- ifelse(data_wide$all_469==1, 1, 
                                  ifelse(data_wide$all_469==2, 2,
                                         ifelse(data_wide$all_675V ==1, 3,
                                                ifelse(data_wide$all_675V ==2, 4, data_wide$all_WT_469))))
data_wide$k13_genotype2 <- ifelse(data_wide$id=="PAT-053", 1, data_wide$k13_genotype2)
data_wide$k13_genotype2 <- ifelse(data_wide$id=="PAT-050", 3, data_wide$k13_genotype2)
data_wide$time_to_assay <- as.numeric(data_wide$time_to_assay)

###################################################
###                    Analyses                 ###
###################################################

#######
# Demographics
######

samples <- subset(data_wide, include==1 | !is.na(LM))

samples %>% 
  group_by(site3) %>%
  summarize(N = n(), median_age = median(Age_yrs, na.rm = TRUE), 
            age_q25 = quantile(Age_yrs,probs = .25,  na.rm = TRUE),
            age_q75 = quantile(Age_yrs,probs = .75,  na.rm = TRUE),
            median_parasit = median(Parasitemia, na.rm = TRUE),
            parasit_q25 = quantile(Parasitemia,probs = .25,  na.rm = TRUE),
            parasit_q75 = quantile(Parasitemia,probs = .75,  na.rm = TRUE),
            time = median(time_to_assay, na.rm = TRUE),
            time_q25 = quantile(time_to_assay, probs = .25, na.rm = TRUE),
            time_q75 = quantile(time_to_assay, probs = .75, na.rm = TRUE),
            coi = median(median_coi, na.rm=TRUE),
            coi_q25 = quantile(median_coi, probs = .25, na.rm = TRUE),
            coi_q75 = quantile(median_coi, probs = .75, na.rm = TRUE))

samples %>% 
  group_by(site3) %>%  
  summarize(coi_mean = median(mean_coi, na.rm=TRUE),
          coi_mean_q25 = quantile(mean_coi, probs = .25, na.rm = TRUE),
          coi_mean_q75 = quantile(mean_coi, probs = .75, na.rm = TRUE))

table(samples$Sex, samples$site3)


wilcox.test(samples$Age_yrs ~ samples$site3)
wilcox.test(samples$Parasitemia ~ samples$site3)
wilcox.test(samples$mean ~ samples$site3)
wilcox.test(samples$time_to_assay ~ samples$site3)
t.test(samples$time_to_assay ~ samples$site3)

ggplot(samples, aes(site3, time_to_assay, fill = site3)) +
  geom_point() + 
  geom_boxplot(alpha =.2, outlier.shape = NA)


chisq.test(samples$Sex, samples$site3)



#######
# Prevalences
######

data_wide_inc <- subset(data_wide, include==1 | !is.na(LM))
table(data_wide_inc$site3)


table(data_wide_inc$all_469Y, data_wide_inc$site3)
fisher.test(data_wide_inc$all_469Y, data_wide_inc$site3)

table(data_wide_inc$all_675V, data_wide_inc$site3)
fisher.test(data_wide_inc$all_675V, data_wide_inc$site3)

table(data_wide_inc$lum_cat, data_wide_inc$site3)
fisher.test(data_wide_inc$lum_cat, data_wide_inc$site3)

table(data_wide_inc$k13.Ala578Ser, data_wide_inc$site3)
fisher.test(data_wide_inc$k13.Ala578Ser, data_wide_inc$site3)

table(data_wide_inc$k13.Gly533Ala, data_wide_inc$site3)
fisher.test(data_wide_inc$k13.Gly533Ala, data_wide_inc$site3)




# which k13 polymorphisms were seen?
k13 <- subset(data_wide_inc, select = str_detect(string = names(data_wide_inc), pattern = "k13"))
snp_there <- colSums(k13, na.rm = TRUE) > 0
names(snp_there)[which(snp_there == T)]
rm(k13, snp_there)


#######
# standard markers
######
st_marker <- subset(data_wide_inc, select = c("id", "site3", "crt.Lys76Thr", "mdr1.Asn86Tyr", "mdr1.Tyr184Phe", "mdr1.Asp1246Tyr", "dhfr.ts.Asn51Ile", "dhfr.ts.Cys59Arg", "dhfr.ts.Ser108Asn", "dhfr.ts.Ile164Leu", "dhps.Ala437Gly", "dhps.Lys540Glu", "dhps.Ala581Gly", "dhps.Ala613Ser", "mdr_CNV", "pm1_CNV", "pm2_CNV", "pm3_CNV"))

st_marker[,3:14][ st_marker[,3:14] == 2 ] <- 1

table(st_marker$crt.Lys76Thr, st_marker$site3)
table(st_marker$mdr1.Asn86Tyr, st_marker$site3)
table(st_marker$mdr1.Tyr184Phe, st_marker$site3)
fisher.test(st_marker$mdr1.Tyr184Phe, st_marker$site3)
table(st_marker$mdr1.Asp1246Tyr, st_marker$site3)
fisher.test(st_marker$mdr1.Asp1246Tyr, st_marker$site3)
table(st_marker$mdr_CNV, st_marker$site3)
table(st_marker$dhfr.ts.Asn51Ile, st_marker$site3)
table(st_marker$dhfr.ts.Cys59Arg, st_marker$site3)
fisher.test(st_marker$dhfr.ts.Cys59Arg, st_marker$site3)
table(st_marker$dhfr.ts.Ser108Asn, st_marker$site3)
table(st_marker$dhfr.ts.Ile164Leu, st_marker$site3)
fisher.test(st_marker$dhfr.ts.Ile164Leu, st_marker$site3)
table(st_marker$dhps.Ala437Gly, st_marker$site3)
fisher.test(st_marker$dhps.Ala437Gly, st_marker$site3)
table(st_marker$dhps.Lys540Glu, st_marker$site3)
fisher.test(st_marker$dhps.Lys540Glu, st_marker$site3)
table(st_marker$dhps.Ala581Gly, st_marker$site3)
fisher.test(st_marker$dhps.Ala581Gly, st_marker$site3)
table(st_marker$pm2_CNV, st_marker$site3)

#######
# Drug sensitivities
######
### Drug sensitivity by region ###
## LM ##
subset(data_wide, !is.na(LM)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))

wilcox.test(LM ~ site3, data_wide, exact = FALSE)

## CQ ##
subset(data_wide, !is.na(CQ)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(CQ, na.rm = TRUE), min = min(CQ, na.rm = TRUE), max = max(CQ, na.rm = TRUE))

wilcox.test(CQ ~ site3, data_wide, exact = FALSE)

## MDAQ ##
subset(data_wide, !is.na(MDAQ)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(MDAQ, na.rm = TRUE), min = min(MDAQ, na.rm = TRUE), max = max(MDAQ, na.rm = TRUE))

wilcox.test(MDAQ ~ site3, data_wide, exact = FALSE)

## PQ ##
subset(data_wide, !is.na(PQ)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(PQ, na.rm = TRUE), min = min(PQ, na.rm = TRUE), max = max(PQ, na.rm = TRUE))

wilcox.test(PQ ~ site3, data_wide, exact = FALSE)

## MQ ##
subset(data_wide, !is.na(MQ)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(MQ, na.rm = TRUE), min = min(MQ, na.rm = TRUE), max = max(MQ, na.rm = TRUE))

wilcox.test(MQ ~ site3, data_wide, exact = FALSE)

## DHA ##
subset(data_wide, !is.na(DHA)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(DHA, na.rm = TRUE), min = min(DHA, na.rm = TRUE), max = max(DHA, na.rm = TRUE))

wilcox.test(DHA ~ site3, data_wide, exact = FALSE)

## PND ##
subset(data_wide, !is.na(PND)) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(PND, na.rm = TRUE), min = min(PND, na.rm = TRUE), max = max(PND, na.rm = TRUE))

wilcox.test(PND ~ site3, data_wide, exact = FALSE)

## RSA ##
subset(data_wide, !is.na(weighted_RSA) & include==1) %>%
  group_by(site3) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))

wilcox.test(weighted_RSA ~ site3, subset(data_wide, include ==1), conf.int = TRUE)

###########
##########
## K13 v Lum ##
subset(data_wide, !is.na(LM)) %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))

wilcox.test(LM ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(LM ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

#data_wide %>%
#  group_by(k13_genotype2) %>%
#  summarize(N = n(), median = median(LM, na.rm = TRUE), min = min(LM, na.rm = TRUE), max = max(LM, na.rm = TRUE))

#wilcox.test(LM ~ k13_genotype2, subset(data_wide, k13_genotype2 == 0 | k13_genotype2==1), conf.int = TRUE)
#wilcox.test(LM ~ k13_genotype2, subset(data_wide, k13_genotype2 == 0 | k13_genotype2==2), conf.int = TRUE)
#wilcox.test(LM ~ k13_genotype2, subset(data_wide, k13_genotype2 == 0 | k13_genotype2==3), conf.int = TRUE)
#wilcox.test(LM ~ k13_genotype2, subset(data_wide, k13_genotype2 == 0 | k13_genotype2==4), conf.int = TRUE)

###########
##########
## K13 v MEF ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(MQ, na.rm = TRUE), min = min(MQ, na.rm = TRUE), max = max(MQ, na.rm = TRUE))

wilcox.test(MQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(MQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

###########
##########
## K13 v CQ ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(CQ, na.rm = TRUE), min = min(CQ, na.rm = TRUE), max = max(CQ, na.rm = TRUE))

wilcox.test(CQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(CQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

###########
##########
## K13 v DHA ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(DHA, na.rm = TRUE), min = min(DHA, na.rm = TRUE), max = max(DHA, na.rm = TRUE))

wilcox.test(DHA ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(DHA ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

###########
##########
## K13 v MDAQ ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(MDAQ, na.rm = TRUE), min = min(MDAQ, na.rm = TRUE), max = max(MDAQ, na.rm = TRUE))

wilcox.test(MDAQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(MDAQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)
###########
##########
## K13 v PQ ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(PQ, na.rm = TRUE), min = min(PQ, na.rm = TRUE), max = max(PQ, na.rm = TRUE))

wilcox.test(PQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(PQ ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)
###########
##########
## K13 v PND ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(PND, na.rm = TRUE), min = min(PND, na.rm = TRUE), max = max(PND, na.rm = TRUE))

wilcox.test(PND ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(PND ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)

###########
##########
## K13 v PYR ##
data_wide %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(PYR, na.rm = TRUE), min = min(PYR, na.rm = TRUE), max = max(PYR, na.rm = TRUE))

wilcox.test(PYR ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==1), conf.int = TRUE)
wilcox.test(PYR ~ k13_genotype, subset(data_wide, k13_genotype == 0 | k13_genotype==2), conf.int = TRUE)
###########
##########
## K13 v RSA ##
# WT includes parasites w/ other K13 mutations

subset(data_wide, include ==1) %>%
  group_by(all_469Y) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))

wilcox.test(weighted_RSA ~ all_469Y, subset(data_wide, include ==1 & (all_469Y == 0 | all_469Y==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ all_469Y, subset(data_wide, include ==1 & (all_469Y == 0 | all_469Y==2)), conf.int = TRUE)

subset(data_wide, include ==1) %>%
  group_by(all_675V) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ all_675V, subset(data_wide, include ==1 & (all_675V == 0 | all_675V==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ all_675V, subset(data_wide, include ==1 & (all_675V == 0 | all_675V==2)), conf.int = TRUE)

# WT vs. Mix/Mut; WT does not include parasites w/ other K13 mutations
subset(data_wide, include ==1) %>%
  group_by(all_WT_469) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ all_WT_469, subset(data_wide, include ==1), conf.int = TRUE)

subset(data_wide, include ==1) %>%
  group_by(all_WT3_469) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ all_WT3_469, subset(data_wide, include ==1 & (all_WT3_469 == 0 | all_WT3_469==1)), conf.int = TRUE, exact = FALSE)
wilcox.test(weighted_RSA ~ all_WT3_469, subset(data_wide, include ==1 & (all_WT3_469 == 0 | all_WT3_469==2)), conf.int = TRUE, exact = FALSE)

subset(data_wide, include ==1) %>%
  group_by(all_WT3_675) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ all_WT3_675, subset(data_wide, include ==1 & (all_WT3_675 == 0 | all_WT3_675==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ all_WT3_675, subset(data_wide, include ==1 & (all_WT3_675 == 0 | all_WT3_675==2)), conf.int = TRUE)
# WT vs. Mix and WT vs. Mut; WT does not include parasites w/ other K13 mutations
subset(data_wide, include ==1) %>%
  group_by(k13_genotype) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))
wilcox.test(weighted_RSA ~ k13_genotype, subset(data_wide, include == 1 & (k13_genotype == 0 | k13_genotype==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype, subset(data_wide, include == 1 & (k13_genotype == 0 | k13_genotype==2)), conf.int = TRUE)


# mix 469, mut 469, mix 675, mut 675 vs wt
subset(data_wide, include ==1) %>%
  group_by(k13_genotype2) %>%
  summarize(N = n(), median = median(weighted_RSA, na.rm = TRUE), min = min(weighted_RSA, na.rm = TRUE), max = max(weighted_RSA, na.rm = TRUE))

wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==1)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==2)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==3)), conf.int = TRUE)
wilcox.test(weighted_RSA ~ k13_genotype2, subset(data_wide, include==1 & (k13_genotype2 == 0 | k13_genotype2==4)), conf.int = TRUE)


data_wide$k13_469_675 <- 
  ifelse(!is.na(data_wide$weighted_RSA) & data_wide$include==1 & data_wide$all_469Y==0 & data_wide$all_675V==0, 0,
  ifelse(!is.na(data_wide$weighted_RSA) & data_wide$include==1 & (data_wide$all_469Y==1 | data_wide$all_469Y==2) & data_wide$all_675V==0, 1,
  ifelse(!is.na(data_wide$weighted_RSA) & data_wide$include==1 & data_wide$all_469Y==0 & (data_wide$all_675V==1 | data_wide$all_675V==2), 2,
  ifelse((!is.na(data_wide$weighted_RSA) & data_wide$include==1 & data_wide$all_469Y==1 | data_wide$all_469Y==2) & (data_wide$all_675V==1 | data_wide$all_675V==2), 3, NA))))

table(data_wide$k13_469_675)       

subset(data_wide, include ==1) %>%
  group_by(k13_469_675) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))

wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 0 | k13_469_675==1)))
wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 0 | k13_469_675==2)))
wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 0 | k13_469_675==3)))
wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 1 | k13_469_675==2)))
wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 1 | k13_469_675==3)))
wilcox.test(weighted_RSA ~ k13_469_675, subset(data_wide, include==1 & (k13_469_675 == 2 | k13_469_675==3)))

