###################################################
###          Install required packages          ###
###################################################
library(sf)
library(raster)
library(dplyr)
library(tmap)
library(leaflet)
library(ggplot2)
library(tidyverse)
library(cartography)
library(ggrepel)

###################################################
###                 Import Data                 ###
###################################################
rm(list=ls())
setwd("C:/Users/conradm1/Desktop/2022 Tumwebaze et al Nature comm")

data(World)
uganda_adm2 <- read_sf("databases/uga_admbnda_ubos_20200824_shp/uga_admbnda_adm2_ubos_20200824.shp")
important_sites <- read.csv("databases/important_sites.csv")
data <- read.csv("databases/PRX-04_k13_data.csv")

###################################################
###                 Tabulate                    ###
###################################################
tab1 <- as.data.frame(xtabs(~anymut + site, data))
tab_wide1 <- spread(tab1, anymut, Freq)
names(tab_wide1) <- c("site", "WT_N", "Mut_N")
tab_wide1$N <- tab_wide1$WT_N + tab_wide1$Mut_N
tab_wide1$Mut_f <- round(tab_wide1$Mut_N/tab_wide1$N * 100, 1)
tab_wide1 <- subset(tab_wide1, select = c("site", "N", "WT_N", "Mut_N", "Mut_f"))

# make categorical value for prevalence of mutations
tab_wide1$cat <- ifelse(tab_wide1$Mut_f==0, 0,
                        ifelse(tab_wide1$Mut_f>0 & tab_wide1$Mut_f<= 10, 1,
                               ifelse(tab_wide1$Mut_f > 10 & tab_wide1$Mut_f <= 20, 2,
                                      ifelse(tab_wide1$Mut_f > 20 & tab_wide1$Mut_f <= 30, 3,
                                             ifelse(tab_wide1$Mut_f > 30 & tab_wide1$Mut_f <= 40, 4, NA)))))

# Merge 2019 Prevalence data w/ map data#
tab_wide1$ADM2_EN <- recode(tab_wide1$site, "AG" = "Agago", "AM" = "Amolatar", "AR" = "Arua", "HO" = "Hoima", "JI" = "Jinja", "KAP" = "Kapchorwa","KB" = "Rukiga", "KBG" = "Kaabong", "KBK" = "Koboko", "KN" = "Kanungu", "KO" = "Kole", "KS" = "Kasese", "KTK" = "Katakwi", "LA" = "Lamwo", "MU" = "Mubende", "TO" = "Tororo")
prx04 <- subset(tab_wide1, select = c("ADM2_EN", "cat"))
prx04_map <- merge(uganda_adm2, prx04, by = "ADM2_EN")

# Format important site data for adding points (sf object)
sites <- st_as_sf(important_sites, coords = c("Longitude", "Latitude"), crs = 4326)

# Subset important sites to only include TDH, Patongo, and Busiu
important_sites_filter <- subset(important_sites, Site=="Tororo" | Site=="Patongo" |Site=="Busiu" | Site=="Masafu")
sites_filter <- subset(sites,  Site=="Tororo" | Site=="Patongo" | Site=="Busiu" | Site=="Masafu")


V1 <-
  ggplot(data = uganda_adm2) +
  geom_sf(fill = 'white', color = "darkgray") +
  geom_sf(data = prx04_map, mapping = aes(fill = as.factor(cat), color = as.factor(cat))) +
  geom_point(data = important_sites_filter, aes(x = Longitude, y = Latitude)) +
  scale_fill_manual(values = c("bisque", "plum", "purple", "darkmagenta","purple4"), labels = c("0%", "1-10%", "11-20%", "20-30%", "30-40%"), name = "Prevalence 469Y & 675V") +  
  scale_color_manual(values = c("bisque", "plum", "purple", "darkmagenta","purple4"), labels = c("0%", "1-10%", "11-20%", "20-30%", "30-40%"), name = "Prevalence 469Y & 675V") +  
  geom_label_repel(data = important_sites_filter, aes(x = Longitude, y = Latitude,label = Site),
                        force = 20, nudge_x = .8, seed = 5) +
  geom_label(aes(label = "Uganda", x = 32.45, y = 1), size = 10) +
  theme_void() +
    theme(legend.position = c(0.7, 0.14), legend.background = element_rect(fill = "white", color = "white", size = 3))
  V1
ggsave(plot = V1, "figures/Exvivo_sites.tiff", width = 5, height = 6)



