###################################################
###          Install required packages          ###
###################################################
library(dplyr)
library(tidyr)
library(stringr)
###################################################
###                 Import Data                 ###
###################################################
#setwd("") #recommend making directory specifically for RSA analysis if IC50 analysis is also planned
# Make new directory within working directory named "Wilcox_output" within working directory where files generated as part of this script will be written

rm(list=ls())
all_data <- read.csv("databases/artR_genes_wide.csv")

#################################################
# subset to include only isolates w/ LUM IC50s
all_rsa <- subset(all_data, !is.na(weighted_RSA) & include==1 & (Month==6 | Month==7) & Year==21)

##########################################
all_rsa_genotypes <- all_rsa[18:1927]
########################
# filter samples
########################
# count missing genotypes
all_rsa$missing_persample <- apply(all_rsa_genotypes, 1, function(x){sum(is.na(x))})
summary(all_rsa$missing_persample)
data_filt_samples <- subset(all_rsa, missing_persample < 1900, select = -c(missing_persample))

########################
# filter loci
########################
filt_rsa_genotypes <- data_filt_samples[18:1927]

# select polymorphic loci (not just WT)
polymorphic_rsa <- colSums(filt_rsa_genotypes, na.rm = TRUE) > 0
table(polymorphic_rsa)

variable_l <- filt_rsa_genotypes[polymorphic_rsa]

# select polymorphic loci (not just Mutant)
polymorphic_rsa_mut <- sapply(variable_l, function(x) any(x != 2, na.rm=T))
table(polymorphic_rsa_mut)
variable_ll <- variable_l[polymorphic_rsa_mut]

variable_lc <- variable_ll
# remove loci that are missing too much data
#missing_perlocus <- sapply(variable_ll, function(x) sum(!is.na(x))) < 69/2
#table(missing_perlocus)
#variable_lc <- variable_ll[!missing_perlocus]

RSA<- subset(data_filt_samples, select= c(id, weighted_RSA, all_WT_469, all_WT_675))
RSA2<- subset(data_filt_samples, select= c(id, weighted_RSA))
RSA3 <- subset(data_filt_samples, select = c(id, weighted_RSA, all_469Y, all_675V))
data_filt <- cbind(RSA, variable_lc)

##########################################
data_long <- pivot_longer(data_filt, Pfubp1.Ala3099Val:PF3D7_1115300.Val147Phe, names_to = "locus", values_to = "genotype")

tab <- as.data.frame(xtabs(~locus + genotype, data_long))
tab_wide <- spread(tab, genotype, Freq)
names(tab_wide) <- c("locus", "WT_N", "Mix_N", "Mut_N")
tab_wide$N <- tab_wide$WT_N + tab_wide$Mix_N + tab_wide$Mut_N
tab_wide$WT_f <- round(tab_wide$WT_N/tab_wide$N * 100, 1)
tab_wide$Mix_f <- round(tab_wide$Mix_N/tab_wide$N * 100, 1)
tab_wide$Mut_f <- round(tab_wide$Mut_N/tab_wide$N * 100, 1)
tab_wide$WT_Mix_f <- round((tab_wide$WT_N + tab_wide$Mix_N)/tab_wide$N * 100, 1)
tab_wide$Mix_Mut_f <- round((tab_wide$Mix_N + tab_wide$Mut_N)/tab_wide$N * 100, 1)

tab_wide <- subset(tab_wide, select = c("locus", "N", "WT_N", "Mix_N", "Mut_N", "WT_f", "Mix_f", "Mut_f",  "WT_Mix_f",  "Mix_Mut_f"))

# higher prevalence of Mut than WT
tab_wide_swap_maj_min <- subset(tab_wide, WT_Mix_f < Mix_Mut_f)
tab_wide_swap_maj_min_high <- subset(tab_wide_swap_maj_min, WT_Mix_f > 15)
tab_wide_swap_maj_min_low <- subset(tab_wide_swap_maj_min, WT_Mix_f <= 15)

loc_swap_maj_min_high <- tab_wide_swap_maj_min_high$locus
loc_swap_maj_min_low <- tab_wide_swap_maj_min_low$locus
swap_maj_min_high <- subset(data_long, locus %in% loc_swap_maj_min_high)
swap_maj_min_low <- subset(data_long, locus %in% loc_swap_maj_min_low)

# For swap, 0 is Mutant (major allele) and 1 is WT or Mixed (minor allele)
swap_maj_min_high$genotype_bi <- ifelse(swap_maj_min_high$genotype==2, 0, 
                                        ifelse(swap_maj_min_high$genotype==0, 1, swap_maj_min_high$genotype))
table(swap_maj_min_high$genotype_bi)
table(swap_maj_min_high$genotype)

swap_maj_min_low$genotype_bi <- ifelse(swap_maj_min_low$genotype==1, 0, swap_maj_min_low$genotype)
swap_maj_min_low$genotype_bi <- ifelse(swap_maj_min_low$genotype_bi==2, 1, swap_maj_min_low$genotype_bi)
table(swap_maj_min_low$genotype_bi)
table(swap_maj_min_low$genotype)

# high prevalence (>15%)
tab_wide_hi_prev <- subset(tab_wide, (WT_Mix_f >= Mix_Mut_f) & Mut_f > 15)
loc_hi_prev <- tab_wide_hi_prev$locus
hi_prev <- subset(data_long, locus %in% loc_hi_prev)
hi_prev$genotype_bi <- ifelse(hi_prev$genotype==2, 1, hi_prev$genotype)

# low prevalence (<=15)
tab_wide_lo_prev <- subset(tab_wide, (WT_Mix_f >= Mix_Mut_f) & Mut_f <= 15)
loc_lo_prev <- tab_wide_lo_prev$locus
lo_prev <- subset(data_long, locus %in% loc_lo_prev)
lo_prev$genotype_bi <- ifelse(lo_prev$genotype==0, 1, lo_prev$genotype)
lo_prev$genotype_bi <- ifelse(lo_prev$genotype_bi==2, 0, lo_prev$genotype_bi)

#rm(list= ls()[!(ls() %in% c('swap_maj_min_high', 'swap_maj_min_low', 'hi_prev', 'lo_prev','RSA','RSA2'))])

#####################################
### make k13_genotype variable
#####################################
# 0 = major genotype, wt 469
# 1 = minor genotype, wt 469
# 2 = major, mutant 469
# 3 = minor genotype, mutant 469
# 4 = major, mutant 675
# 5 = minor genotype, mutant 675
hi_prev$k13_genotype <- 
  ifelse(hi_prev$genotype_bi==1 & hi_prev$all_WT_469==1, 3,
   ifelse(hi_prev$genotype_bi==0 & hi_prev$all_WT_469==1, 2,
    ifelse(hi_prev$genotype_bi==1  & hi_prev$all_WT_469==0, 1,
     ifelse(hi_prev$all_WT_469==0  & hi_prev$all_WT_469==0, 0, NA))))
hi_prev$k13_genotype <- 
  ifelse(is.na(hi_prev$k13_genotype) & hi_prev$genotype_bi==0 & hi_prev$all_WT_675==1, 4,
    ifelse(is.na(hi_prev$k13_genotype) & hi_prev$genotype_bi==1 & hi_prev$all_WT_675==1, 5, hi_prev$k13_genotype))
table(hi_prev$k13_genotype)

swap_maj_min_high$k13_genotype <- 
  ifelse(swap_maj_min_high$genotype_bi==1 & swap_maj_min_high$all_WT_469==1, 3,
    ifelse(swap_maj_min_high$genotype_bi== 0 & swap_maj_min_high$all_WT_469==1, 2,
      ifelse(swap_maj_min_high$genotype_bi==1  & swap_maj_min_high$all_WT_469==0, 1,
        ifelse(swap_maj_min_high$all_WT_469==0 & swap_maj_min_high$all_WT_469==0, 0, NA))))
swap_maj_min_high$k13_genotype <- 
  ifelse(is.na(swap_maj_min_high$k13_genotype) & swap_maj_min_high$genotype_bi==0 & swap_maj_min_high$all_WT_675==1, 4,
         ifelse(is.na(swap_maj_min_high$k13_genotype) & swap_maj_min_high$genotype_bi==1 & swap_maj_min_high$all_WT_675==1, 5, swap_maj_min_high$k13_genotype))
table(swap_maj_min_high$k13_genotype)


#hi_prev 
#swap_maj_min_high
#################
# N and medians
#################
# WT, Mix, Mut
hi_prev_wmm <- hi_prev %>% 
  group_by(locus, genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(hi_prev_wmm, "Wilcox_output/output_medians_hi_prev_wmm.csv", row.names = F)

swap_hi_wmm <- swap_maj_min_high %>% 
  group_by(locus, genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(swap_hi_wmm, "Wilcox_output/output_medians_swap_hi_wmm.csv", row.names = F)

# Minority present vs minority absent
hi_prev_wm <- hi_prev %>% 
  group_by(locus, genotype_bi) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(hi_prev_wm, "Wilcox_output/output_medians_hi_prev_wm.csv", row.names = F)

swap_hi_wm <- swap_maj_min_high %>% 
  group_by(locus, genotype_bi) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(swap_hi_wm, "Wilcox_output/output_medians_swap_hi_wm.csv", row.names = F)

# W/ k13 mutations
hi_prev_k13wm <- hi_prev %>% 
  group_by(locus, k13_genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(hi_prev_k13wm, "Wilcox_output/output_medians_hi_prev_k13wm.csv", row.names = F)

swap_hi_k13wm <- swap_maj_min_high %>% 
  group_by(locus, k13_genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(swap_hi_k13wm, "Wilcox_output/output_medians_swap_hi_k13wm.csv", row.names = F)

######################################
#  Wilcox Minority vs No Minority
######################################
# High prevalence
hi_prev_wm_wide <- subset(hi_prev, select = c("id", "locus", "genotype_bi"))
hi_prev_wm_wide <- pivot_wider(hi_prev_wm_wide, names_from = locus, values_from = genotype_bi)
hi_prev_wm_wide2 <- merge(RSA2, hi_prev_wm_wide, by = "id")

column_names <- data.frame(locus = "locus", 
                           N_samples = "N", 
                           N_minority = "N_minority", 
                           N_no_minority = "N_no_minority", 
                           median_minority = "median_minority", 
                           median_no_minority = "median_no_minority", 
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_hi_prev_wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:32){ #j = loci
  locus<-colnames(hi_prev_wm_wide2)[j]
  N_samples <- length(which((hi_prev_wm_wide2[,j]==0 | hi_prev_wm_wide2[,j]==1)))
  N_minority <- length(which((hi_prev_wm_wide2[,j]==1)))
  N_no_minority <- length(which((hi_prev_wm_wide2[,j]==0)))
  median_minority <- median(hi_prev_wm_wide2$weighted_RSA[which(hi_prev_wm_wide2[j]==1)], na.rm=TRUE)
  median_no_minority <- median(hi_prev_wm_wide2$weighted_RSA[which(hi_prev_wm_wide2[j]==0)], na.rm=TRUE)
  tes <- wilcox.test(hi_prev_wm_wide2$weighted_RSA ~ hi_prev_wm_wide2[,j])$p.value
  newline <- data.frame(t(c(locus, N_samples, N_minority, N_no_minority, median_minority, median_no_minority, tes)))
  write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}
##########
# High prevalence where minority allele is the WT allele
swap_hi_wm_wide <- subset(swap_maj_min_high, select = c("id", "locus", "genotype_bi"))
swap_hi_wm_wide <- pivot_wider(swap_hi_wm_wide, names_from = locus, values_from = genotype_bi)
swap_hi_wm_wide2 <- merge(RSA2, swap_hi_wm_wide, by = "id")

column_names <- data.frame(locus = "locus", 
                           N_samples = "N", 
                           N_minority = "N_minority", 
                           N_no_minority = "N_no_minority", 
                           median_minority = "median_minority", 
                           median_no_minority = "median_no_minority", 
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_swap_hi_wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:24){ 
  locus<-colnames(swap_hi_wm_wide2)[j]
  N_samples <- length(which((swap_hi_wm_wide2[,j]==0 | swap_hi_wm_wide2[,j]==1)))
  N_minority <- length(which((swap_hi_wm_wide2[,j]==1)))
  N_no_minority <- length(which((swap_hi_wm_wide2[,j]==0)))
  median_minority <- median(swap_hi_wm_wide2$weighted_RSA[which(swap_hi_wm_wide2[j]==-1)], na.rm=TRUE)
  median_no_minority <- median(swap_hi_wm_wide2$weighted_RSA[which(swap_hi_wm_wide2[j]==0)], na.rm=TRUE)
  tes <- wilcox.test(swap_hi_wm_wide2$weighted_RSA ~ swap_hi_wm_wide2[,j])$p.value
  newline <- data.frame(t(c(locus, N_samples, N_minority, N_no_minority, median_minority, median_no_minority, tes)))
  write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}

#########################################
# WT vs Mix, WT vs. Mut
##########################################
# High prevalence
hi_prev_wmm_wide <- subset(hi_prev, select = c("id", "locus", "genotype"))
hi_prev_wmm_wide <- pivot_wider(hi_prev_wmm_wide, names_from = locus, values_from = genotype)
hi_prev_wmm_wide2 <- merge(RSA2, hi_prev_wmm_wide, by = "id")

column_names2 <- data.frame(locus = "locus", 
                            N_samples = "N", 
                            N_wt = "N_wt", 
                            N_mix = "N_mix",
                            N_mut = "N_mut",
                            medianwt = "median_wt", 
                            medianmix = "median_mix",
                            medianmut = "median_mut",
                            tes1 = "MutvMix_p-value",
                            tes2 = "MutvWt_p-value")
write.table(column_names2, file = "Wilcox_output/output_wilcox_hi_prev_wmm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:32){ 
  if(any(hi_prev_wmm_wide2[,j]==0, na.rm=T) & any(hi_prev_wmm_wide2[,j]==1, na.rm=T)  & any(hi_prev_wmm_wide2[,j]==2, na.rm=T)) 
  {
    locus<-colnames(hi_prev_wmm_wide2)[j]
    N_samples <- length(which((hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1 | hi_prev_wmm_wide2[,j]==2)))
    N_wt <- length(which((hi_prev_wmm_wide2[,j]==0)))
    N_mix <-  length(which((hi_prev_wmm_wide2[,j]==1)))
    N_mut <-  length(which((hi_prev_wmm_wide2[,j]==2)))
    medianwt <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==0)], na.rm=TRUE)
    medianmix <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==1)], na.rm=TRUE)
    medianmut <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==2)], na.rm=TRUE)
    medianmut
    temp1 <- hi_prev_wmm_wide2[which(hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1),]
    temp2 <- hi_prev_wmm_wide2[which(hi_prev_wmm_wide2[,j]== 0| hi_prev_wmm_wide2[,j]==2),]
    tes1 <- wilcox.test(temp1$weighted_RSA ~ temp1[,j])$p.value
    tes2 <- wilcox.test(temp2$weighted_RSA ~ temp2[,j])$p.value
    newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    if(any(hi_prev_wmm_wide2[,j]==0, na.rm=T) & any(hi_prev_wmm_wide2[,j]==1, na.rm=T)) 
    {
      locus<-colnames(hi_prev_wmm_wide2)[j]
      N_samples <- length(which((hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1 | hi_prev_wmm_wide2[,j]==2)))
      N_wt <- length(which((hi_prev_wmm_wide2[,j]==0)))
      N_mix <-  length(which((hi_prev_wmm_wide2[,j]==1)))
      N_mut <-  length(which((hi_prev_wmm_wide2[,j]==2)))
      medianwt <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==0)], na.rm=TRUE)
      medianmix <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==1)], na.rm=TRUE)
      medianmut <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==2)], na.rm=TRUE)
      temp1 <- hi_prev_wmm_wide2[which(hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1),]
      tes1 <- wilcox.test(temp1$weighted_RSA ~ temp1[,j])$p.value
      tes2 <- "NA"
      newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
      write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    }
    else {
      if (any(hi_prev_wmm_wide2[,j]==0, na.rm=T) & any(hi_prev_wmm_wide2[,j]==2, na.rm=T)) 
      {
        locus<-colnames(hi_prev_wmm_wide2)[j]
        N_samples <- length(which((hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1 | hi_prev_wmm_wide2[,j]==2)))
        N_wt <- length(which((hi_prev_wmm_wide2[,j]==0)))
        N_mix <-  length(which((hi_prev_wmm_wide2[,j]==1)))
        N_mut <-  length(which((hi_prev_wmm_wide2[,j]==2)))
        medianwt <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==0)], na.rm=TRUE)
        medianmix <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==1)], na.rm=TRUE)
        medianmut <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==2)], na.rm=TRUE)
        temp2 <- hi_prev_wmm_wide2[which(hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==2),]
        tes2 <- wilcox.test(temp2$weighted_RSA ~ temp2[,j])$p.value
        tes1 <- "NA"
        newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
        write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
      }
      else {
        locus <-colnames(hi_prev_wmm_wide2)[j]
        N_samples <- length(which((hi_prev_wmm_wide2[,j]==0 | hi_prev_wmm_wide2[,j]==1 | hi_prev_wmm_wide2[,j]==2)))
        N_wt <- length(which((hi_prev_wmm_wide2[,j]==0)))
        N_mix <-  length(which((hi_prev_wmm_wide2[,j]==1)))
        N_mut <-  length(which((hi_prev_wmm_wide2[,j]==2)))
        medianwt <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==0)], na.rm=TRUE)
        medianmix <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==1)], na.rm=TRUE)
        medianmut <- median(hi_prev_wmm_wide2$weighted_RSA[which(hi_prev_wmm_wide2[j]==2)], na.rm=TRUE)
        tes1 <- "NA"
        tes2 <- "NA"
        newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
        write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
      }
    }
  }
}

############################################################
# High prevalence where minority allele is the WT allele
############################################################
swap_hi_wmm_wide <- subset(swap_maj_min_high, select = c("id", "locus", "genotype"))
swap_hi_wmm_wide <- pivot_wider(swap_hi_wmm_wide, names_from = locus, values_from = genotype)
swap_hi_wmm_wide2 <- merge(RSA2, swap_hi_wmm_wide, by = "id")

column_names2 <- data.frame(locus = "locus", 
                            N_samples = "N", 
                            N_wt = "N_wt", 
                            N_mix = "N_mix",
                            N_mut = "N_mut",
                            medianwt = "median_wt", 
                            medianmix = "median_mix",
                            medianmut = "median_mut",
                            tes1 = "MutvMix_p-value",
                            tes2 = "MutvWT_p-value")
write.table(column_names2, file = "Wilcox_output/output_wilcox_swap_hi_wmm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:24){ 
  if(any(swap_hi_wmm_wide2[,j]==0, na.rm=T) & any(swap_hi_wmm_wide2[,j]==1, na.rm=T)  & any(swap_hi_wmm_wide2[,j]==2, na.rm=T)) 
  {
    locus<-colnames(swap_hi_wmm_wide2)[j]
    N_samples <- length(which((swap_hi_wmm_wide2[,j]==0 | swap_hi_wmm_wide2[,j]==1 | swap_hi_wmm_wide2[,j]==2)))
    N_wt <- length(which((swap_hi_wmm_wide2[,j]==0)))
    N_mix <-  length(which((swap_hi_wmm_wide2[,j]==1)))
    N_mut <-  length(which((swap_hi_wmm_wide2[,j]==2)))
    medianwt <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==0)], na.rm=TRUE)
    medianmix <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==1)], na.rm=TRUE)
    medianmut <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==2)], na.rm=TRUE)
    temp1 <- swap_hi_wmm_wide2[which(swap_hi_wmm_wide2[,j]==2 | swap_hi_wmm_wide2[,j]==1),]
    temp2 <- swap_hi_wmm_wide2[which(swap_hi_wmm_wide2[,j]==2 | swap_hi_wmm_wide2[,j]==0),]
    tes1 <- wilcox.test(temp1$weighted_RSA ~ temp1[,j])$p.value
    tes2 <- wilcox.test(temp2$weighted_RSA ~ temp2[,j])$p.value
    newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    if(any(swap_hi_wmm_wide2[,j]==0, na.rm=T) & any(swap_hi_wmm_wide2[,j]==1, na.rm=T)) 
    {
      locus<-colnames(swap_hi_wmm_wide2)[j]
      N_samples <- length(which((swap_hi_wmm_wide2[,j]==0 | swap_hi_wmm_wide2[,j]==1 | swap_hi_wmm_wide2[,j]==2)))
      N_wt <- length(which((swap_hi_wmm_wide2[,j]==0)))
      N_mix <-  length(which((swap_hi_wmm_wide2[,j]==1)))
      N_mut <-  length(which((swap_hi_wmm_wide2[,j]==2)))
      medianwt <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==0)], na.rm=TRUE)
      medianmix <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==1)], na.rm=TRUE)
      medianmut <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==2)], na.rm=TRUE)
      temp1 <- swap_hi_wmm_wide2[which(swap_hi_wmm_wide2[,j]==2 | swap_hi_wmm_wide2[,j]==1),]
      tes1 <- wilcox.test(temp1$weighted_RSA ~ temp1[,j])$p.value
      tes2 <- ""
      newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
      write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    }
    else {
      if (any(swap_hi_wmm_wide2[,j]==0, na.rm=T) & any(swap_hi_wmm_wide2[,j]==2, na.rm=T)) 
      {
        locus<-colnames(swap_hi_wmm_wide2)[j]
        N_samples <- length(which((swap_hi_wmm_wide2[,j]==0 | swap_hi_wmm_wide2[,j]==1 | swap_hi_wmm_wide2[,j]==2)))
        N_wt <- length(which((swap_hi_wmm_wide2[,j]==0)))
        N_mix <-  length(which((swap_hi_wmm_wide2[,j]==1)))
        N_mut <-  length(which((swap_hi_wmm_wide2[,j]==2)))
        medianwt <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==0)], na.rm=TRUE)
        medianmix <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==1)], na.rm=TRUE)
        medianmut <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==2)], na.rm=TRUE)
        temp2 <- swap_hi_wmm_wide2[which(swap_hi_wmm_wide2[,j]==2 | swap_hi_wmm_wide2[,j]==0),]
        tes2 <- wilcox.test(temp2$weighted_RSA ~ temp2[,j])$p.value
        tes1 <- ""
        newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
        write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
      }
      else {
        locus <-colnames(swap_hi_wmm_wide2)[j]
        N_samples <- length(which((swap_hi_wmm_wide2[,j]==0 | swap_hi_wmm_wide2[,j]==1 | swap_hi_wmm_wide2[,j]==2)))
        N_wt <- length(which((swap_hi_wmm_wide2[,j]==0)))
        N_mix <-  length(which((swap_hi_wmm_wide2[,j]==1)))
        N_mut <-  length(which((swap_hi_wmm_wide2[,j]==2)))
        medianwt <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==0)], na.rm=TRUE)
        medianmix <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==1)], na.rm=TRUE)
        medianmut <- median(swap_hi_wmm_wide2$weighted_RSA[which(swap_hi_wmm_wide2[j]==2)], na.rm=TRUE)
        tes1 <- ""
        tes2 <- ""
        newline <- data.frame(t(c(locus, N_samples, N_wt, N_mix, N_mut, medianwt, medianmix, medianmut, tes1, tes2)))
        write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_wmm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
      }
    }
  }
}

##########################################
# K13 vs mutants
##########################################
# High prevalence
hi_prev_k13wm_wide <- subset(hi_prev, select = c("id", "locus", "k13_genotype"))
hi_prev_k13wm_wide <- pivot_wider(hi_prev_k13wm_wide, names_from = locus, values_from = k13_genotype)
hi_prev_k13wm_wide2 <- merge(RSA2, hi_prev_k13wm_wide, by = "id")

column_names <- data.frame(locus = "locus", 
                           comparison = "comparison",
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:32){   
  if(any(hi_prev_k13wm_wide2[,j]==0, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==1, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "1_K13wt-wt vs. K13wt-mut (0 v 1)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==0 | hi_prev_k13wm_wide2[,j]==1),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "2_K13wt-wt vs. K13wt-mut (0 v 1)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==0, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==2, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "3_K13wt-wt vs. 469Y-wt (0 v 2)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==0 | hi_prev_k13wm_wide2[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "4_K13wt-wt vs. 469Y-wt (0 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==0, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "5_K13wt-wt vs. 469Y-mut (0 v 3)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==0 | hi_prev_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "6_K13wt-wt vs. 469Y-mut (0 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==0, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==4, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "7_K13wt-wt vs. 675V-wt (0 v 4)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==0 | hi_prev_k13wm_wide2[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "8_K13wt-wt vs. 675V-wt (0 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==0, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "9_K13wt-wt vs. 675V-mut (0 v 5)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==0 | hi_prev_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "10_K13wt-wt vs. 675V-mut (0 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  
  if(any(hi_prev_k13wm_wide2[,j]==1, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==2, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <- "11_K13wt-mut vs. 469Y-wt (1 v 2)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==1 | hi_prev_k13wm_wide2[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "12_K13wt-mut vs. 469Y-wt (1 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==1, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "13_K13wt-mut vs. 469Y-mut (1 v 3)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==1 | hi_prev_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "14_K13wt-mut vs. 469Y-mut (1 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==2, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "15_469Y-wt vs. 469Y-mut (2 v 3)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==2 | hi_prev_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "16_469Y-wt vs. 469Y-mut (2 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==1, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==4, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "17_K13wt-mut vs. 675V-wt (1 v 4)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==1 | hi_prev_k13wm_wide2[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "18_K13wt-mut vs. 675V-wt (1 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==1, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <- "19_K13wt-mut vs. 675V-mut (1 v 5)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==1 | hi_prev_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "20_K13wt-mut vs. 675V-mut (1 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(hi_prev_k13wm_wide2[,j]==4, na.rm=T) & any(hi_prev_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(hi_prev_k13wm_wide2)[j]
    comparison <-  "21_675V-wt vs. 675V-mut (4 v 5)"
    temp <- hi_prev_k13wm_wide2[which(hi_prev_k13wm_wide2[,j]==4 | hi_prev_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(hi_prev_wmm_wide2)[j]
    comparison <-  "22_675V-wt vs. 675V-mut (4 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_hi_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
}

# High prevalence where minority allele is the WT allele
swap_hi_k13wm_wide <- subset(swap_maj_min_high, select = c("id", "locus", "k13_genotype"))
swap_hi_k13wm_wide <- pivot_wider(swap_hi_k13wm_wide, names_from = locus, values_from = k13_genotype)
swap_hi_k13wm_wide2 <- merge(RSA2, swap_hi_k13wm_wide, by = "id")

column_names <- data.frame(locus = "locus", 
                           comparison = "comparison",
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:24){   
  if(any(swap_hi_k13wm_wide2[,j]==0, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==1, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "1_K13wt-wt vs. K13wt-mut (0 v 1)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==0 | swap_hi_k13wm_wide2[,j]==1),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "2_K13wt-wt vs. K13wt-mut (0 v 1)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==0, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==2, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "3_K13wt-wt vs. 469Y-wt (0 v 2)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==0 | swap_hi_k13wm_wide2[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "4_K13wt-wt vs. 469Y-wt (0 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==0, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "5_K13wt-wt vs. 469Y-mut (0 v 3)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==0 | swap_hi_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "6_K13wt-wt vs. 469Y-mut (0 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==0, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==4, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "7_K13wt-wt vs. 675V-wt (0 v 4)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==0 | swap_hi_k13wm_wide2[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "8_K13wt-wt vs. 675V-wt (0 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==0, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "9_K13wt-wt vs. 675V-mut (0 v 5)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==0 | swap_hi_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "10_K13wt-wt vs. 675V-mut (0 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  
  if(any(swap_hi_k13wm_wide2[,j]==1, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==2, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <- "11_K13wt-mut vs. 469Y-wt (1 v 2)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==1 | swap_hi_k13wm_wide2[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "12_K13wt-mut vs. 469Y-wt (1 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==1, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "13_K13wt-mut vs. 469Y-mut (1 v 3)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==1 | swap_hi_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "14_K13wt-mut vs. 469Y-mut (1 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==2, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==3, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "15_469Y-wt vs. 469Y-mut (2 v 3)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==2 | swap_hi_k13wm_wide2[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "16_469Y-wt vs. 469Y-mut (2 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==1, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==4, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "17_K13wt-mut vs. 675V-wt (1 v 4)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==1 | swap_hi_k13wm_wide2[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "18_K13wt-mut vs. 675V-wt (1 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==1, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <- "19_K13wt-mut vs. 675V-mut (1 v 5)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==1 | swap_hi_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "20_K13wt-mut vs. 675V-mut (1 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(swap_hi_k13wm_wide2[,j]==4, na.rm=T) & any(swap_hi_k13wm_wide2[,j]==5, na.rm=T)) 
  {
    locus <-colnames(swap_hi_k13wm_wide2)[j]
    comparison <-  "21_675V-wt vs. 675V-mut (4 v 5)"
    temp <- swap_hi_k13wm_wide2[which(swap_hi_k13wm_wide2[,j]==4 | swap_hi_k13wm_wide2[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(swap_hi_wmm_wide2)[j]
    comparison <-  "22_675V-wt vs. 675V-mut (4 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_swap_hi_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
}
##########################################
###               lo prev: WT vs Mut
##########################################
#rm(list= ls()[!(ls() %in% c('lo_prev', 'swap_maj_min_low', 'RSA', 'RSA2'))])

# subset by gene
acs10_1 <- subset(lo_prev, str_detect(locus, "0525100"), select = c(id, locus, genotype))
acs10_2 <- subset(swap_maj_min_low, str_detect(locus, "0525100"), select = c(id, locus, genotype))
ap2mu <- subset(lo_prev, str_detect(locus, "Pfap2mu"), select = c(id, locus, genotype))
ap3d_1 <- subset(lo_prev, str_detect(locus, "0808100"), select = c(id, locus, genotype))
ap3d_2 <- subset(swap_maj_min_low, str_detect(locus, "0808100"), select = c(id, locus, genotype))
arps10 <- subset(lo_prev, str_detect(locus, "1460900"), select = c(id, locus, genotype))
atp6 <- subset(lo_prev, str_detect(locus, "atp6"), select = c(id, locus, genotype))
coronin <- subset(lo_prev,  str_detect(locus, "1251200"), select = c(id, locus, genotype))
crt <- subset(lo_prev, str_detect(locus, "crt"), select = c(id, locus, genotype))
fd <- subset(lo_prev, str_detect(locus, "fd"), select = c(id, locus, genotype))
fp2a <- subset(lo_prev, str_detect(locus, "1115700"), select = c(id, locus, genotype))
fp2b <- subset(lo_prev, str_detect(locus, "1115300"), select = c(id, locus, genotype))
fp3 <- subset(lo_prev, str_detect(locus, "1115400"), select = c(id, locus, genotype))
k13 <- subset(lo_prev, str_detect(locus, "k13"), select = c(id, locus, genotype))
kelch10 <- subset(lo_prev, str_detect(locus, "kelch10"), select = c(id, locus, genotype))
mcp <- subset(lo_prev, str_detect(locus, "mcp"), select = c(id, locus, genotype))
mdr2 <- subset(lo_prev, str_detect(locus, "mdr2"), select = c(id, locus, genotype))
mrp1 <- subset(lo_prev, str_detect(locus, "mrp1"), select = c(id, locus, genotype))
mrp2_1 <- subset(lo_prev, str_detect(locus, "pfmrp2"), select = c(id, locus, genotype))
mrp2_2 <- subset(swap_maj_min_low, str_detect(locus, "pfmrp2"), select = c(id, locus, genotype))
PF3D7_0806800 <- subset(lo_prev, str_detect(locus, "0806800"), select = c(id, locus, genotype))
PF3D7_0808300 <- subset(lo_prev, str_detect(locus, "0808300"), select = c(id, locus, genotype))
PF3D7_1322700 <- subset(lo_prev, str_detect(locus, "1322700"), select = c(id, locus, genotype))
PF3D7_1433800 <- subset(lo_prev, str_detect(locus, "1433800"), select = c(id, locus, genotype))
PF3D7_1451200 <- subset(lo_prev, str_detect(locus, "1451200"), select = c(id, locus, genotype))
PI4K <- subset(lo_prev, str_detect(locus, "PI4K"), select = c(id, locus, genotype))
pib7 <- subset(lo_prev, str_detect(locus, "pib7"), select = c(id, locus, genotype))
pph <- subset(lo_prev, str_detect(locus, "pph"), select = c(id, locus, genotype))
rad14 <- subset(lo_prev, str_detect(locus, "0710400"), select = c(id, locus, genotype))
rpn10 <- subset(lo_prev, str_detect(locus, "0807800"), select = c(id, locus, genotype))
ruvb2 <- subset(lo_prev, str_detect(locus, "1106000"), select = c(id, locus, genotype))
sec14 <- subset(lo_prev, str_detect(locus, "0626400"), select = c(id, locus, genotype))
ubp <- subset(lo_prev, str_detect(locus, "Pfubp1"), select = c(id, locus, genotype))
ubp2 <- subset(swap_maj_min_low, str_detect(locus, "Pfubp1"), select = c(id, locus, genotype))

acs10 <- rbind(acs10_1, acs10_2)
ap3d <- rbind(ap3d_1, ap3d_2)
mrp2 <- rbind(mrp2_1, mrp2_2)

# convert to wide format
acs10_wide <- pivot_wider(data = acs10, names_from = locus, values_from = genotype)
ap2mu_wide <- pivot_wider(data = ap2mu, names_from = locus, values_from = genotype)
ap3d_wide <- pivot_wider(data = ap3d, names_from = locus, values_from = genotype)
arps10_wide <- pivot_wider(data = arps10, names_from = locus, values_from = genotype)
atp6_wide <- pivot_wider(data = atp6, names_from = locus, values_from = genotype)
coronin_wide <- pivot_wider(data = coronin, names_from = locus, values_from = genotype)
crt_wide <- pivot_wider(data = crt, names_from = locus, values_from = genotype)
fd_wide <- pivot_wider(data = fd, names_from = locus, values_from = genotype)
fp2a_wide <- pivot_wider(data = fp2a, names_from = locus, values_from = genotype)
fp2b_wide <- pivot_wider(data = fp2b, names_from = locus, values_from = genotype)
fp3_wide <- pivot_wider(data = fp3, names_from = locus, values_from = genotype)
k13_wide <- pivot_wider(data = k13, names_from = locus, values_from = genotype)
kelch10_wide <- pivot_wider(data = kelch10, names_from = locus, values_from = genotype)
mcp_wide <- pivot_wider(data = mcp, names_from = locus, values_from = genotype)
mdr2_wide <- pivot_wider(data = mdr2, names_from = locus, values_from = genotype)
mrp1_wide <- pivot_wider(data = mrp1, names_from = locus, values_from = genotype)
mrp2_wide <- pivot_wider(data = mrp2, names_from = locus, values_from = genotype)
PF3D7_0806800_wide <- pivot_wider(data = PF3D7_0806800, names_from = locus, values_from = genotype)
PF3D7_0808300_wide <- pivot_wider(data = PF3D7_0808300, names_from = locus, values_from = genotype)
PF3D7_1322700_wide <- pivot_wider(data = PF3D7_1322700, names_from = locus, values_from = genotype)
PF3D7_1433800_wide <- pivot_wider(data = PF3D7_1433800, names_from = locus, values_from = genotype)
PF3D7_1451200_wide <- pivot_wider(data = PF3D7_1451200, names_from = locus, values_from = genotype)
PI4K_wide <- pivot_wider(data = PI4K, names_from = locus, values_from = genotype)
pib7_wide <- pivot_wider(data = pib7, names_from = locus, values_from = genotype)
pph_wide <- pivot_wider(data = pph, names_from = locus, values_from = genotype)
rad14_wide <- pivot_wider(data = rad14, names_from = locus, values_from = genotype)
rpn10_wide <- pivot_wider(data = rpn10, names_from = locus, values_from = genotype)
ruvb2_wide <- pivot_wider(data = ruvb2, names_from = locus, values_from = genotype)
sec14_wide <- pivot_wider(data = sec14, names_from = locus, values_from = genotype)
ubp_wide <- pivot_wider(data = ubp, names_from = locus, values_from = genotype)

# determine number of variables
length_acs10 <- length(names(acs10_wide))
length_ap2mu <- length(names(ap2mu_wide))
length_ap3d <- length(names(ap3d_wide))
length_arps10 <- length(names(arps10_wide))
length_atp6 <- length(names(atp6_wide))
length_coronin <- length(names(coronin_wide))
length_crt <- length(names(crt_wide))
length_fd <- length(names(fd_wide))
length_fp2a <- length(names(fp2a_wide))
length_fp2b <- length(names(fp2b_wide))
length_fp3 <- length(names(fp3_wide))
length_k13 <- length(names(k13_wide))
length_kelch10 <- length(names(kelch10_wide))
length_mcp <- length(names(mcp_wide))
length_mdr2 <- length(names(mdr2_wide))
length_mrp1 <- length(names(mrp1_wide))
length_mrp2 <- length(names(mrp2_wide))
length_PF3D7_0806800 <- length(names(PF3D7_0806800_wide))
length_PF3D7_0808300 <- length(names(PF3D7_0808300_wide))
length_PF3D7_1322700 <- length(names(PF3D7_1322700_wide))
length_PF3D7_1433800 <- length(names(PF3D7_1433800_wide))
length_PF3D7_1451200 <- length(names(PF3D7_1451200_wide))
length_PI4K <- length(names(PI4K_wide))
length_pib7 <- length(names(pib7_wide))
length_pph <- length(names(pph_wide))
length_rad14 <- length(names(rad14_wide))
length_rpn10 <- length(names(rpn10_wide))
length_ruvb2 <- length(names(ruvb2_wide))
length_sec14 <- length(names(sec14_wide))
length_ubp <- length(names(ubp_wide))

# names of loci for each gene
loci_acs10 <- names(acs10_wide)
loci_ap2mu <- names(ap2mu_wide)
loci_ap3d <- names(ap3d_wide)
loci_arps10 <- names(arps10_wide)
loci_atp6 <- names(atp6_wide)
loci_coronin <- names(coronin_wide)
loci_crt <- names(crt_wide)
loci_fd <- names(fd_wide)
loci_fp2a <- names(fp2a_wide)
loci_fp2b <- names(fp2b_wide)
loci_fp3 <- names(fp3_wide)
loci_k13 <- names(k13_wide)
loci_kelch10 <- names(kelch10_wide)
loci_mcp <- names(mcp_wide)
loci_mdr2 <- names(mdr2_wide)
loci_mrp1 <- names(mrp1_wide)
loci_mrp2 <- names(mrp2_wide)
loci_PF3D7_0806800 <- names(PF3D7_0806800_wide)
loci_PF3D7_0808300 <- names(PF3D7_0808300_wide)
loci_PF3D7_1322700 <- names(PF3D7_1322700_wide)
loci_PF3D7_1433800 <- names(PF3D7_1433800_wide)
loci_PF3D7_1451200 <- names(PF3D7_1451200_wide)
loci_PI4K <- names(PI4K_wide)
loci_pib7 <- names(pib7_wide)
loci_pph <- names(pph_wide)
loci_rad14 <- names(rad14_wide)
loci_rpn10 <- names(rpn10_wide)
loci_ruvb2 <- names(ruvb2_wide)
loci_sec14 <- names(sec14_wide)
loci_ubp <- names(ubp_wide)


all_loci <- c(loci_acs10,loci_ap2mu, loci_ap3d, loci_arps10, loci_atp6, loci_coronin, loci_crt, loci_fd, loci_fp2a, loci_fp2b, loci_fp3, loci_k13, loci_kelch10, loci_mcp, loci_mdr2, loci_mrp1, loci_mrp2, loci_PF3D7_0806800, loci_PF3D7_0808300, loci_PF3D7_1322700, loci_PF3D7_1433800, loci_PF3D7_1451200, loci_PI4K, loci_pib7, loci_pph, loci_rad14, loci_rpn10, loci_ruvb2, loci_sec14, loci_ubp)
write.csv(all_loci, "Wilcox_output/merged_files/all_loci.csv")

# tabulate loci with only 1 variable site
fp2_rsa <- merge(fp3_wide, RSA2, by = "id")
table(fp2_rsa$PF3D7_1115400.Asn468Tyr)
fp2_rsa %>% 
  group_by(PF3D7_1115400.Asn468Tyr) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(fp2_rsa$weighted_RSA[which(fp2_rsa$PF3D7_1115400.Asn468Tyr==0 | fp2_rsa$PF3D7_1115400.Asn468Tyr==1)] ~ fp2_rsa$PF3D7_1115400.Asn468Tyr[which(fp2_rsa$PF3D7_1115400.Asn468Tyr==0 | fp2_rsa$PF3D7_1115400.Asn468Tyr==1)])
wilcox.test(fp2_rsa$weighted_RSA[which(fp2_rsa$PF3D7_1115400.Asn468Tyr==0 | fp2_rsa$PF3D7_1115400.Asn468Tyr==2)] ~ fp2_rsa$PF3D7_1115400.Asn468Tyr[which(fp2_rsa$PF3D7_1115400.Asn468Tyr==0 | fp2_rsa$PF3D7_1115400.Asn468Tyr==2)])


kelch10_rsa <- merge(kelch10_wide, RSA2, by = "id")
table(kelch10_rsa$kelch10.Leu622Ile)
kelch10_rsa %>% 
  group_by(kelch10.Leu622Ile) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(kelch10_rsa$weighted_RSA[which(kelch10_rsa$kelch10.Leu622Ile==0 | kelch10_rsa$kelch10.Leu622Ile==1)] ~ kelch10_rsa$kelch10.Leu622Ile[which(kelch10_rsa$kelch10.Leu622Ile==0 |kelch10_rsa$kelch10.Leu622Ile==1)])


mcp_rsa <- merge(mcp_wide, RSA2, by = "id")
table(mcp_rsa$mcp.Pro280Thr)
mcp_rsa %>% 
  group_by(mcp.Pro280Thr) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(mcp_rsa$weighted_RSA[which(mcp_rsa$mcp.Pro280Thr==0 | mcp_rsa$mcp.Pro280Thr==1)] ~ mcp_rsa$mcp.Pro280Thr[which(mcp_rsa$mcp.Pro280Thr==0 | mcp_rsa$mcp.Pro280Thr==1)])
wilcox.test(mcp_rsa$weighted_RSA[which(mcp_rsa$mcp.Pro280Thr==0 | mcp_rsa$mcp.Pro280Thr==2)] ~ mcp_rsa$mcp.Pro280Thr[which(mcp_rsa$mcp.Pro280Thr==0 | mcp_rsa$mcp.Pro280Thr==2)])


mdr2_rsa <- merge(mdr2_wide, RSA2, by = "id")
table(mdr2_rsa$mdr2.Ile492Val)
mdr2_rsa %>% 
  group_by(mdr2.Ile492Val) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(mdr2_rsa$weighted_RSA[which(mdr2_rsa$mdr2.Ile492Val==0 | mdr2_rsa$mdr2.Ile492Val==1)] ~ mdr2_rsa$mdr2.Ile492Val[which(mdr2_rsa$mdr2.Ile492Val==0 | mdr2_rsa$mdr2.Ile492Val==1)])
wilcox.test(mdr2_rsa$weighted_RSA[which(mdr2_rsa$mdr2.Ile492Val==0 | mdr2_rsa$mdr2.Ile492Val==2)] ~ mdr2_rsa$mdr2.Ile492Val[which(mdr2_rsa$mdr2.Ile492Val==0 | mdr2_rsa$mdr2.Ile492Val==2)])


PI4K_rsa <- merge(PI4K_wide, RSA2, by = "id")
table(PI4K_rsa$PI4K.His886Leu)
PI4K_rsa %>% 
  group_by(PI4K.His886Leu) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(PI4K_rsa$weighted_RSA[which(PI4K_rsa$PI4K.His886Leu==0 | PI4K_rsa$PI4K.His886Leu==1)] ~ PI4K_rsa$PI4K.His886Leu[which(PI4K_rsa$PI4K.His886Leu==0 | PI4K_rsa$PI4K.His886Leu==1)])

PF3D7_0808300_rsa <- merge(PF3D7_0808300_wide, RSA2, by = "id")
table(PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile)
PF3D7_0808300_rsa %>% 
  group_by(PF3D7_0808300.Met144Ile) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(PF3D7_0808300_rsa$weighted_RSA[which(PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==0 | PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==1)] ~ PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile[which(PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==0 | PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==1)])
wilcox.test(PF3D7_0808300_rsa$weighted_RSA[which(PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==0 | PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==2)] ~ PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile[which(PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==0 | PF3D7_0808300_rsa$PF3D7_0808300.Met144Ile==2)])

# make "any" variable
acs10_wide$any_acs10 <- ifelse(rowSums(acs10_wide[,c(2:length_acs10)], na.rm = TRUE) > 0 &  rowSums(!is.na(acs10_wide[,c(2:length_acs10)])) > 0, 1,  ifelse(rowSums(acs10_wide[,c(2:length_acs10)], na.rm = TRUE) == 0 &   rowSums(!is.na(acs10_wide[,c(2:length_acs10)])) > 0, 0, NA))

ap2mu_wide$any_ap2mu <- ifelse(rowSums(ap2mu_wide[,c(2:length_ap2mu)], na.rm = TRUE) > 0 &  rowSums(!is.na(ap2mu_wide[,c(2:length_ap2mu)])) > 0, 1,  ifelse(rowSums(ap2mu_wide[,c(2:length_ap2mu)], na.rm = TRUE) == 0 &   rowSums(!is.na(ap2mu_wide[,c(2:length_ap2mu)])) > 0, 0, NA))

ap3d_wide$any_ap3d <- ifelse(rowSums(ap3d_wide[,c(2:length_ap3d)], na.rm = TRUE) > 0 &  rowSums(!is.na(ap3d_wide[,c(2:length_ap3d)])) > 0, 1,  ifelse(rowSums(ap3d_wide[,c(2:length_ap3d)], na.rm = TRUE) == 0 &   rowSums(!is.na(ap3d_wide[,c(2:length_ap3d)])) > 0, 0, NA))

#arps10_wide$any_arps10 <- ifelse(rowSums(arps10_wide[,c(2:length_arps10)], na.rm = TRUE) > 0 & rowSums(!is.na(arps10_wide[,c(2:length_arps10)])) > 0, 1, ifelse(rowSums(arps10_wide[,c(2:length_arps10)], na.rm = TRUE) == 0 & rowSums(!is.na(arps10_wide[,c(2:length_arps10)])) > 0, 0, NA))

atp6_wide$any_atp6 <- ifelse(rowSums(atp6_wide[,c(2:length_atp6)], na.rm = TRUE) > 0 &  rowSums(!is.na(atp6_wide[,c(2:length_atp6)])) > 0, 1,  ifelse(rowSums(atp6_wide[,c(2:length_atp6)], na.rm = TRUE) == 0 &   rowSums(!is.na(atp6_wide[,c(2:length_atp6)])) > 0, 0, NA))

coronin_wide$any_coronin <- ifelse(rowSums(coronin_wide[,c(2:length_coronin)], na.rm = TRUE) > 0 &  rowSums(!is.na(coronin_wide[,c(2:length_coronin)])) > 0, 1,  ifelse(rowSums(coronin_wide[,c(2:length_coronin)], na.rm = TRUE) == 0 &   rowSums(!is.na(coronin_wide[,c(2:length_coronin)])) > 0, 0, NA))

crt_wide$any_crt <- ifelse(rowSums(crt_wide[,c(2:length_crt)], na.rm = TRUE) > 0 &  rowSums(!is.na(crt_wide[,c(2:length_crt)])) > 0, 1,  ifelse(rowSums(crt_wide[,c(2:length_crt)], na.rm = TRUE) == 0 &   rowSums(!is.na(crt_wide[,c(2:length_crt)])) > 0, 0, NA))

#fd_wide$any_fd <- ifelse(rowSums(fd_wide[,c(2:length_fd)], na.rm = TRUE) > 0 &  rowSums(!is.na(fd_wide[,c(2:length_fd)])) > 0, 1,  ifelse(rowSums(fd_wide[,c(2:length_fd)], na.rm = TRUE) == 0 &   rowSums(!is.na(fd_wide[,c(2:length_fd)])) > 0, 0, NA))

fp2a_wide$any_fp2a <- ifelse(rowSums(fp2a_wide[,c(2:length_fp2a)], na.rm = TRUE) > 0 &  rowSums(!is.na(fp2a_wide[,c(2:length_fp2a)])) > 0, 1,  ifelse(rowSums(fp2a_wide[,c(2:length_fp2a)], na.rm = TRUE) == 0 &   rowSums(!is.na(fp2a_wide[,c(2:length_fp2a)])) > 0, 0, NA))

fp2b_wide$any_fp2b <- ifelse(rowSums(fp2b_wide[,c(2:length_fp2b)], na.rm = TRUE) > 0 &  rowSums(!is.na(fp2b_wide[,c(2:length_fp2b)])) > 0, 1,  ifelse(rowSums(fp2b_wide[,c(2:length_fp2b)], na.rm = TRUE) == 0 &   rowSums(!is.na(fp2b_wide[,c(2:length_fp2b)])) > 0, 0, NA))

fp3_wide$any_fp3 <- ifelse(rowSums(fp3_wide[,c(2:length_fp3)], na.rm = TRUE) > 0 &  rowSums(!is.na(fp3_wide[,c(2:length_fp3)])) > 0, 1,  ifelse(rowSums(fp3_wide[,c(2:length_fp3)], na.rm = TRUE) == 0 &   rowSums(!is.na(fp3_wide[,c(2:length_fp3)])) > 0, 0, NA))

k13_wide$any_k13 <- ifelse(rowSums(k13_wide[,c(2:length_k13)], na.rm = TRUE) > 0 &  rowSums(!is.na(k13_wide[,c(2:length_k13)]))> 0, 1,  ifelse(rowSums(k13_wide[,c(2:length_k13)], na.rm = TRUE) == 0 &   rowSums(!is.na(k13_wide[,c(2:length_k13)])) > 0, 0, NA))

kelch10_wide$any_kelch10 <- ifelse(rowSums(kelch10_wide[,c(2:length_kelch10)], na.rm = TRUE) > 0 &  rowSums(!is.na(kelch10_wide[,c(2:length_kelch10)])) > 0, 1,  ifelse(rowSums(kelch10_wide[,c(2:length_kelch10)], na.rm = TRUE) == 0 &   rowSums(!is.na(kelch10_wide[,c(2:length_kelch10)])) > 0, 0, NA))

mcp_wide$any_mcp <- ifelse(rowSums(mcp_wide[,c(2:length_mcp)], na.rm = TRUE) > 0 &  rowSums(!is.na(mcp_wide[,c(2:length_mcp)])) > 0, 1,  ifelse(rowSums(mcp_wide[,c(2:length_mcp)], na.rm = TRUE) == 0 &   rowSums(!is.na(mcp_wide[,c(2:length_mcp)])) > 0, 0, NA))

mdr2_wide$any_mdr2 <- ifelse(rowSums(mdr2_wide[,c(2:length_mdr2)], na.rm = TRUE) > 0 &  rowSums(!is.na(mdr2_wide[,c(2:length_mdr2)])) > 0, 1,  ifelse(rowSums(mdr2_wide[,c(2:length_mdr2)], na.rm = TRUE) == 0 &   rowSums(!is.na(mdr2_wide[,c(2:length_mdr2)])) > 0, 0, NA))

mrp1_wide$any_mrp1 <- ifelse(rowSums(mrp1_wide[,c(2:length_mrp1)], na.rm = TRUE) > 0 &  rowSums(!is.na(mrp1_wide[,c(2:length_mrp1)])) > 0, 1,  ifelse(rowSums(mrp1_wide[,c(2:length_mrp1)], na.rm = TRUE) == 0 &   rowSums(!is.na(mrp1_wide[,c(2:length_mrp1)])) > 0, 0, NA))

mrp2_wide$any_mrp2 <- ifelse(rowSums(mrp2_wide[,c(2:length_mrp2)], na.rm = TRUE) > 0 &  rowSums(!is.na(mrp2_wide[,c(2:length_mrp2)])) > 0, 1,  ifelse(rowSums(mrp2_wide[,c(2:length_mrp2)], na.rm = TRUE) == 0 &   rowSums(!is.na(mrp2_wide[,c(2:length_mrp2)])) > 0, 0, NA))

PF3D7_0806800_wide$any_PF3D7_0806800 <- ifelse(rowSums(PF3D7_0806800_wide[,c(2:length_PF3D7_0806800)], na.rm = TRUE) > 0 &  rowSums(!is.na(PF3D7_0806800_wide[,c(2:length_PF3D7_0806800)])) > 0, 1,  ifelse(rowSums(PF3D7_0806800_wide[,c(2:length_PF3D7_0806800)], na.rm = TRUE) == 0 &   rowSums(!is.na(PF3D7_0806800_wide[,c(2:length_PF3D7_0806800)])) > 0, 0, NA))

PF3D7_0808300_wide$any_PF3D7_0808300 <- ifelse(rowSums(PF3D7_0808300_wide[,c(2:length_PF3D7_0808300)], na.rm = TRUE) > 0 &  rowSums(!is.na(PF3D7_0808300_wide[,c(2:length_PF3D7_0808300)])) > 0, 1,  ifelse(rowSums(PF3D7_0808300_wide[,c(2:length_PF3D7_0808300)], na.rm = TRUE) == 0 &   rowSums(!is.na(PF3D7_0808300_wide[,c(2:length_PF3D7_0808300)])) > 0, 0, NA))

PF3D7_1322700_wide$any_PF3D7_1322700 <- ifelse(rowSums(PF3D7_1322700_wide[,c(2:length_PF3D7_1322700)], na.rm = TRUE) > 0 &  rowSums(!is.na(PF3D7_1322700_wide[,c(2:length_PF3D7_1322700)])) > 0, 1,  ifelse(rowSums(PF3D7_1322700_wide[,c(2:length_PF3D7_1322700)], na.rm = TRUE) == 0 &   rowSums(!is.na(PF3D7_1322700_wide[,c(2:length_PF3D7_1322700)])) > 0, 0, NA))

PF3D7_1433800_wide$any_PF3D7_1433800 <- ifelse(rowSums(PF3D7_1433800_wide[,c(2:length_PF3D7_1433800)], na.rm = TRUE) > 0 &  rowSums(!is.na(PF3D7_1433800_wide[,c(2:length_PF3D7_1433800)])) > 0, 1,  ifelse(rowSums(PF3D7_1433800_wide[,c(2:length_PF3D7_1433800)], na.rm = TRUE) == 0 &   rowSums(!is.na(PF3D7_1433800_wide[,c(2:length_PF3D7_1433800)])) > 0, 0, NA))

PF3D7_1451200_wide$any_PF3D7_1451200 <- ifelse(rowSums(PF3D7_1451200_wide[,c(2:length_PF3D7_1451200)], na.rm = TRUE) > 0 &  rowSums(!is.na(PF3D7_1451200_wide[,c(2:length_PF3D7_1451200)])) > 0, 1,  ifelse(rowSums(PF3D7_1451200_wide[,c(2:length_PF3D7_1451200)], na.rm = TRUE) == 0 &   rowSums(!is.na(PF3D7_1451200_wide[,c(2:length_PF3D7_1451200)])) > 0, 0, NA))

PI4K_wide$any_PI4K <- ifelse(rowSums(PI4K_wide[,c(2:length_PI4K)], na.rm = TRUE) > 0 &  rowSums(!is.na(PI4K_wide[,c(2:length_PI4K)])) > 0, 1,  ifelse(rowSums(PI4K_wide[,c(2:length_PI4K)], na.rm = TRUE) == 0 &   rowSums(!is.na(PI4K_wide[,c(2:length_PI4K)])) > 0, 0, NA))

#pib7_wide$any_pib7 <- ifelse(rowSums(pib7_wide[,c(2:length_pib7)], na.rm = TRUE) > 0 &  rowSums(!is.na(pib7_wide[,c(2:length_pib7)])) > 0, 1,  ifelse(rowSums(pib7_wide[,c(2:length_pib7)], na.rm = TRUE) == 0 &   rowSums(!is.na(pib7_wide[,c(2:length_pib7)])) > 0, 0, NA))

pph_wide$any_pph <- ifelse(rowSums(pph_wide[,c(2:length_pph)], na.rm = TRUE) > 0 &  rowSums(!is.na(pph_wide[,c(2:length_pph)])) > 0, 1,  ifelse(rowSums(pph_wide[,c(2:length_pph)], na.rm = TRUE) == 0 &   rowSums(!is.na(pph_wide[,c(2:length_pph)])) > 0, 0, NA))

rad14_wide$any_rad14 <- ifelse(rowSums(rad14_wide[,c(2:length_rad14)], na.rm = TRUE) > 0 &  rowSums(!is.na(rad14_wide[,c(2:length_rad14)])) > 0, 1,  ifelse(rowSums(rad14_wide[,c(2:length_rad14)], na.rm = TRUE) == 0 &   rowSums(!is.na(rad14_wide[,c(2:length_rad14)])) > 0, 0, NA))

rpn10_wide$any_rpn10 <- ifelse(rowSums(rpn10_wide[,c(2:length_rpn10)], na.rm = TRUE) > 0 &  rowSums(!is.na(rpn10_wide[,c(2:length_rpn10)])) > 0, 1,  ifelse(rowSums(rpn10_wide[,c(2:length_rpn10)], na.rm = TRUE) == 0 &   rowSums(!is.na(rpn10_wide[,c(2:length_rpn10)])) > 0, 0, NA))

ruvb2_wide$any_ruvb2 <- ifelse(rowSums(ruvb2_wide[,c(2:length_ruvb2)], na.rm = TRUE) > 0 &  rowSums(!is.na(ruvb2_wide[,c(2:length_ruvb2)])) > 0, 1,  ifelse(rowSums(ruvb2_wide[,c(2:length_ruvb2)], na.rm = TRUE) == 0 &   rowSums(!is.na(ruvb2_wide[,c(2:length_ruvb2)])) > 0, 0, NA))

#sec14_wide$any_sec14 <- ifelse(rowSums(sec14_wide[,c(2:length_sec14)], na.rm = TRUE) > 0 &  rowSums(!is.na(sec14_wide[,c(2:length_sec14)])) > 0, 1,  ifelse(rowSums(sec14_wide[,c(2:length_sec14)], na.rm = TRUE) == 0 &   rowSums(!is.na(sec14_wide[,c(2:length_sec14)])) > 0, 0, NA))

# too polymorphic; need to break it up into segments
#ubp_wide$any_ubp <- ifelse(rowSums(ubp_wide[,c(2:length_ubp)], na.rm = TRUE) > 0 &  rowSums(!is.na(ubp_wide[,c(2:length_ubp)])) > 0, 1,  ifelse(rowSums(ubp_wide[,c(2:length_ubp)], na.rm = TRUE) == 0 &   rowSums(!is.na(ubp_wide[,c(2:length_ubp)])) > 0, 0, NA))

# make subset of "any" variables for merging
acs10_2 <- subset(acs10_wide, select = c(id, any_acs10))
ap2mu_2 <- subset(ap2mu_wide, select = c(id, any_ap2mu))
ap3d_2 <- subset(ap3d_wide, select = c(id, any_ap3d))
#arps10_2 <- subset(arps10_wide, select = c(id, any_arps10))
atp6_2 <- subset(atp6_wide, select = c(id, any_atp6))
coronin_2 <- subset(coronin_wide, select = c(id, any_coronin))
crt_2 <- subset(crt_wide, select = c(id, any_crt))
#fd_2 <- subset(fd_wide, select = c(id, any_fd))
fp2a_2 <- subset(fp2a_wide, select = c(id, any_fp2a))
fp2b_2 <- subset(fp2b_wide, select = c(id, any_fp2b))
fp3_2 <- subset(fp3_wide, select = c(id, any_fp3))
k13_2 <- subset(k13_wide, select = c(id, any_k13))
kelch10_2 <- subset(kelch10_wide, select = c(id, any_kelch10))
mcp_2 <- subset(mcp_wide, select = c(id, any_mcp))
mdr2_2 <- subset(mdr2_wide, select = c(id, any_mdr2))
mrp1_2 <- subset(mrp1_wide, select = c(id, any_mrp1))
mrp2_2 <- subset(mrp2_wide, select = c(id, any_mrp2))
PF3D7_0806800_2 <- subset(PF3D7_0806800_wide, select = c(id, any_PF3D7_0806800))
PF3D7_0808300_2 <- subset(PF3D7_0808300_wide, select = c(id, any_PF3D7_0808300))
PF3D7_1322700_2 <- subset(PF3D7_1322700_wide, select = c(id, any_PF3D7_1322700))
PF3D7_1433800_2 <- subset(PF3D7_1433800_wide, select = c(id, any_PF3D7_1433800))
PF3D7_1451200_2 <- subset(PF3D7_1451200_wide, select = c(id, any_PF3D7_1451200))
PI4K_2 <- subset(PI4K_wide, select = c(id, any_PI4K))
#pib7_2 <- subset(pib7_wide, select = c(id, any_pib7))
pph_2 <- subset(pph_wide, select = c(id, any_pph))
rad14_2 <- subset(rad14_wide, select = c(id, any_rad14))
rpn10_2 <- subset(rpn10_wide, select = c(id, any_rpn10))
ruvb2_2 <- subset(ruvb2_wide, select = c(id, any_ruvb2))
#sec14_2 <- subset(sec14_wide, select = c(id, any_sec14))
#ubp_2 <- subset(ubp_wide, select = c(id, any_ubp))

# merge "any variables"
data_merge <- merge(acs10_2, ap2mu_2, by = "id")
data_merge <- merge(data_merge, ap3d_2, by = "id")
data_merge <- merge(data_merge, atp6_2, by = "id")
data_merge <- merge(data_merge, coronin_2, by = "id")
data_merge <- merge(data_merge, crt_2, by = "id")
data_merge <- merge(data_merge, fp2a_2, by = "id")
data_merge <- merge(data_merge, fp2b_2, by = "id")
data_merge <- merge(data_merge, fp3_2, by = "id")
data_merge <- merge(data_merge, k13_2, by = "id")
data_merge <- merge(data_merge, kelch10_2, by = "id")
data_merge <- merge(data_merge, mcp_2, by = "id")
data_merge <- merge(data_merge, mdr2_2, by = "id")
data_merge <- merge(data_merge, mrp1_2, by = "id")
data_merge <- merge(data_merge, mrp2_2, by = "id")
data_merge <- merge(data_merge, PF3D7_0806800_2, by = "id")
data_merge <- merge(data_merge, PF3D7_0808300_2, by = "id")
data_merge <- merge(data_merge, PF3D7_1322700_2, by = "id")
data_merge <- merge(data_merge, PF3D7_1433800_2, by = "id")
data_merge <- merge(data_merge, PF3D7_1451200_2, by = "id")
data_merge <- merge(data_merge, PI4K_2, by = "id")
data_merge <- merge(data_merge, pph_2, by = "id")
data_merge <- merge(data_merge, rad14_2, by = "id")
data_merge <- merge(data_merge, rpn10_2, by = "id")
data_merge <- merge(data_merge, ruvb2_2, by = "id")

data_merge_bi <- merge(RSA2, data_merge, by = "id")

data_merge_biK13 <- merge(RSA3, data_merge, by = "id")

#######################




#rm(list= ls()[!(ls() %in% c('lo_prev', "data_merge", "data_merge_bi", 'RSA', 'RSA2'))])
###################################
###  Swap to long format
###################################
data_merge_long <- pivot_longer(data_merge, cols = any_acs10:any_ruvb2, names_to = "locus", values_to = "genotype")
data_merge_k13 <- merge(RSA, data_merge_long, by = "id")

data_merge_k13$k13_genotype <- 
  ifelse(data_merge_k13$genotype==1 & data_merge_k13$all_WT_469==1, 3,
         ifelse(data_merge_k13$genotype== 0 & data_merge_k13$all_WT_469==1, 2,
                ifelse(data_merge_k13$genotype==1  & data_merge_k13$all_WT_469==0, 1,
                       ifelse(data_merge_k13$all_WT_469==0 & data_merge_k13$all_WT_469==0, 0, NA))))
data_merge_k13$k13_genotype <- 
  ifelse(is.na(data_merge_k13$k13_genotype) & data_merge_k13$genotype==0 & data_merge_k13$all_WT_675==1, 4,
         ifelse(is.na(data_merge_k13$k13_genotype) & data_merge_k13$genotype==1 & data_merge_k13$all_WT_675==1, 5, data_merge_k13$k13_genotype))
table(data_merge_k13$k13_genotype)

data_merge_k13_ss <- subset(data_merge_k13, select = c(id, locus, k13_genotype))
data_k13geno <- pivot_wider(data_merge_k13_ss, names_from = locus, values_from = k13_genotype)
data_goiK13 <- merge(RSA, data_k13geno, by = "id")
data_goiK13$any_k132 <- ifelse(data_goiK13$all_WT_469==1, NA,
                               ifelse(data_goiK13$all_WT_675==1, NA, data_goiK13$any_k13))

table(data_goiK13$any_k13)
table(data_goiK13$any_k132)

data_goiK13 %>% 
  group_by(any_k132) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
wilcox.test(data_goiK13$weighted_RSA ~ data_goiK13$any_k132)

#################
#  Wilcox WT v Mut
#################

lo_prev_wm <- data_merge_k13 %>% 
  group_by(locus, genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(lo_prev_wm, "Wilcox_output/output_medians_lo_prev_wm.csv", row.names = F)


column_names <- data.frame(locus = "locus",
                           N_samples = "N_samples",
                           N_wt = "N_wt",
                           N_mut = "N_mut",
                           median_wt = "median_wt",
                           median_mut = "median_mut",
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_lo_prev_wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:27){ #j = swapci
  locus <- colnames(data_merge_bi)[j]
  N_samples <- length(which((data_merge_bi[,j]==0 | data_merge_bi[,j]==1)))
  N_wt <- length(which(data_merge_bi[,j]==0))
  N_mut <- length(which(data_merge_bi[,j]==1))
  median_wt <- median(data_merge_bi$weighted_RSA[which(data_merge_bi[j]==0)], na.rm=TRUE)
  median_mut <- median(data_merge_bi$weighted_RSA[which(data_merge_bi[j]==1)], na.rm=TRUE)
  tes<- wilcox.test(data_merge_bi$weighted_RSA ~ data_merge_bi[,j])$p.value
  newline <- data.frame(t(c(locus, N_samples, N_wt, N_mut, median_wt, median_mut, tes)))
  write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
}


##########################################
##########################################
###               data_goiK13
##########################################
##########################################

# 0 = wt/mix genotype, wt 469
# 1 = mutant genotype, wt 469
# 2 = wt/mix genotype, mutant 469
# 3 = mutant genotype, mutant 469
# 4 = wt/mix genotype, mutant 675
# 5 = mutant genotype, mutant 675

lo_prev_k13wm <- data_merge_k13 %>% 
  group_by(locus, k13_genotype) %>%
  summarize(N = n(),
            RSA_median = median(weighted_RSA, na.rm=TRUE),
            RSA_q25 = quantile(weighted_RSA, probs = .25, na.rm = TRUE),
            RSA_q75 = quantile(weighted_RSA, probs = .75, na.rm = TRUE))
write.csv(lo_prev_k13wm, "Wilcox_output/output_medians_lo_prev_k13wm.csv", row.names = F)

##########################################
# K13 vs mutants
##########################################
# 0 = wt genotype, wt 469
# 1 = mutant genotype, wt 469
# 2 = wt genotype, mutant 469
# 3 = mutant genotype, mutant 469
# 4 = wt genotype, mutant 675
# 5 = mutant genotype, mutant 675

column_names <- data.frame(locus = "locus", 
                           comparison = "comparison",
                           tes = "p-value")
write.table(column_names, file = "Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",", quote = TRUE)

for(j in 3:29){   
  if(any(data_goiK13[,j]==0, na.rm=T) & any(data_goiK13[,j]==1, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. K13wt-mut (0 v 1)"
    temp <- data_goiK13[which(data_goiK13[,j]==0 | data_goiK13[,j]==1),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. K13wt-mut (0 v 1)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==0, na.rm=T) & any(data_goiK13[,j]==2, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 469Y-wt (0 v 2)"
    temp <- data_goiK13[which(data_goiK13[,j]==0 | data_goiK13[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 469Y-wt (0 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==0, na.rm=T) & any(data_goiK13[,j]==3, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 469Y-mut (0 v 3)"
    temp <- data_goiK13[which(data_goiK13[,j]==0 | data_goiK13[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 469Y-mut (0 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==0, na.rm=T) & any(data_goiK13[,j]==4, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 675V-wt (0 v 4)"
    temp <- data_goiK13[which(data_goiK13[,j]==0 | data_goiK13[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
    
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 675V-wt (0 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==0, na.rm=T) & any(data_goiK13[,j]==5, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 675V-mut (0 v 5)"
    temp <- data_goiK13[which(data_goiK13[,j]==0 | data_goiK13[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-wt vs. 675V-mut (0 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  
  if(any(data_goiK13[,j]==1, na.rm=T) & any(data_goiK13[,j]==2, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <- "K13wt-mut vs. 469Y-wt (1 v 2)"
    temp <- data_goiK13[which(data_goiK13[,j]==1 | data_goiK13[,j]==2),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 469Y-wt (1 v 2)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==1, na.rm=T) & any(data_goiK13[,j]==3, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 469Y-mut (1 v 3)"
    temp <- data_goiK13[which(data_goiK13[,j]==1 | data_goiK13[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 469Y-mut (1 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==2, na.rm=T) & any(data_goiK13[,j]==3, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "469Y-wt vs. 469Y-mut (2 v 3)"
    temp <- data_goiK13[which(data_goiK13[,j]==2 | data_goiK13[,j]==3),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "469Y-wt vs. 469Y-mut (2 v 3)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==1, na.rm=T) & any(data_goiK13[,j]==4, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 675V-wt (1 v 4)"
    temp <- data_goiK13[which(data_goiK13[,j]==1 | data_goiK13[,j]==4),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 675V-wt (1 v 4)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==1, na.rm=T) & any(data_goiK13[,j]==5, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <- "K13wt-mut vs. 675V-mut (1 v 5)"
    temp <- data_goiK13[which(data_goiK13[,j]==1 | data_goiK13[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "K13wt-mut vs. 675V-mut (1 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  if(any(data_goiK13[,j]==4, na.rm=T) & any(data_goiK13[,j]==5, na.rm=T)) 
  {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "675V-wt vs. 675V-mut (4 v 5)"
    temp <- data_goiK13[which(data_goiK13[,j]==4 | data_goiK13[,j]==5),]
    tes <- (wilcox.test(temp$weighted_RSA ~ temp[,j])$p.value)
    newline <- data.frame(t(c(locus, comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
  else {
    locus <-colnames(data_goiK13)[j]
    comparison <-  "675V-wt vs. 675V-mut (4 v 5)"
    tes <- "NA"
    newline <- data.frame(t(c(locus,comparison, tes)))
    write.table(newline, file="Wilcox_output/output_wilcox_lo_prev_k13wm.csv", row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",")
  }
}



