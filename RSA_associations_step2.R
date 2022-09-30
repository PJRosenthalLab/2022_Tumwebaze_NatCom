# Change working directory to location of output from RSA_associations.R file
#setwd("Wilcox_output/")
#make directory called "merged_files" where output will be written

#rm(list=ls())

## Medians##
hi_med_wm <- read.csv("output_medians_hi_prev_wm.csv")
hi_med_k13wm <- read.csv("output_medians_hi_prev_k13wm.csv")
hi_med_wmm <- read.csv("output_medians_hi_prev_wmm.csv")

swap_hi_med_wm <- read.csv("output_medians_swap_hi_wm.csv")
swap_hi_med_k13wm <- read.csv("output_medians_swap_hi_k13wm.csv")
swap_hi_med_wmm <- read.csv("output_medians_swap_hi_wmm.csv")

lo_med_wm <- read.csv("output_medians_lo_prev_wm.csv")
lo_med_k13wm <- read.csv("output_medians_lo_prev_k13wm.csv")

# Statistics ##
hi_wil_wm <- read.csv("output_wilcox_hi_prev_wm.csv")
hi_wil_k13wm <- read.csv("output_wilcox_hi_prev_k13wm.csv")
hi_wil_wmm <- read.csv("output_wilcox_hi_prev_wmm.csv")

swap_hi_wil_wm <- read.csv("output_wilcox_swap_hi_wm.csv")
swap_hi_wil_k13wm <- read.csv("output_wilcox_swap_hi_k13wm.csv")
swap_hi_wil_wmm<- read.csv("output_wilcox_swap_hi_wmm.csv")

lo_wil_wm <- read.csv("output_wilcox_lo_prev_wm.csv")
lo_wil_k13wm <- read.csv("output_wilcox_lo_prev_k13wm.csv")

#####################
# combine median files
#####################
# Binomial genotypes
hi_med_wm_ss <- subset(hi_med_wm, locus != "k13.Cys469Tyr" & !is.na(genotype_bi))
lo_med_wm_ss <- subset(lo_med_wm, !is.na(genotype))
names(lo_med_wm_ss)[2] <- 'genotype_bi'
swap_med_wm_ss <- subset(swap_hi_med_wm, !is.na(genotype_bi))

medians_wm <- rbind(hi_med_wm_ss, lo_med_wm_ss, swap_med_wm_ss)

# Trinomial genotypes
hi_med_wmm_ss <- subset(hi_med_wmm, locus != "k13.Cys469Tyr" & !is.na(genotype))
swap_med_wmm_ss <- subset(swap_hi_med_wmm, !is.na(genotype))

medians_wmm <- rbind(hi_med_wmm_ss, swap_med_wmm_ss)

# K13 genotypes
hi_med_k13wm_ss <- subset(hi_med_k13wm, locus != "k13.Cys469Tyr" & !is.na(k13_genotype))
lo_med_k13wm_ss <- subset(lo_med_k13wm, !is.na(k13_genotype))
swap_med_k13wm_ss <- subset(swap_hi_med_k13wm, !is.na(k13_genotype))

medians_k13wm <- rbind(hi_med_k13wm_ss, lo_med_k13wm_ss,swap_med_k13wm_ss)

######################
# combine p-values
#######################
names(lo_wil_wm) <- c("locus", "N", "N_no_minority", "N_minority", "median_no_minority", "median_minority", "p.value")
wil_wm <- rbind(hi_wil_wm, lo_wil_wm, swap_hi_wil_wm)
wil_wm_ss <- subset(wil_wm, select = c(locus, N_minority, N_no_minority, p.value))
wil_wm_long <- pivot_longer(wil_wm_ss, cols = N_minority:N_no_minority, names_to = "genotype", values_to = "N")
wil_wm_long$genotype_bi <- ifelse(wil_wm_long$genotype=="N_minority", 1, 
                                  ifelse(wil_wm_long$genotype=="N_no_minority", 0, NA))

# Trinomial genotypes
wil_wmm <- rbind(hi_wil_wmm, swap_hi_wil_wmm)
wil_wmm_ss <- subset(wil_wmm, select = c(locus, N_wt, N_mix, N_mut, WTvMix_p.value, WTvMut_p.value))
wil_wmm_long <- pivot_longer(wil_wmm_ss, cols = N_wt:N_mut, names_to = "genotype_cat", values_to = "N")
wil_wmm_long$genotype <- ifelse(wil_wmm_long$genotype_cat=="N_wt", 0, 
                                  ifelse(wil_wmm_long$genotype_cat=="N_mix", 1,
                                         ifelse(wil_wmm_long$genotype_cat=="N_mut", 2, NA)))
# K13 genotypes
wil_k13wm <- rbind(hi_wil_k13wm, lo_wil_k13wm, swap_hi_wil_k13wm)


##############################
# merge medians and p-values
##############################
wm <- merge(medians_wm, wil_wm_long, by = c("locus", "genotype_bi"), all = T)
wmm <- merge(medians_wmm, wil_wmm_long, by = c("locus", "genotype"), all = T)

write.csv(wm, "merged_files/wm.csv", row.names = F)
write.csv(wmm, "merged_files/wmm.csv", row.names = F)
