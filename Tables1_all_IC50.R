###################################################
###          Install required packages          ###
###################################################
library("tidyverse")
library("lubridate")
library("dplyr")

###################################################
###            Import required files            ###
###################################################
rm(list=ls())

#setwd("")

data_all <- read.csv("databases/exvivo_all_wide.csv")
data <- subset(data_all, Year < 21 | (Year==21 & Month < 8) | (Year==21 & Month==8 & Day <=20))


data %>% 
  summarize(N = n(),
            Age_median = median(Age_yrs, na.rm=TRUE),
            Age_q25 = quantile(Age_yrs, probs = .25, na.rm = TRUE),
            Age_q75 = quantile(Age_yrs, probs = .75, na.rm = TRUE))

table(data$Sex)

data %>% 
  summarize(N = n(),
           Parasitemia_median = median(Parasitemia, na.rm=TRUE),
           Parasitemia_q25 = quantile(Parasitemia, probs = .25, na.rm = TRUE),
           Parasitemia_q75 = quantile(Parasitemia, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            COI_median = median(median, na.rm=TRUE),
            COI_q25 = quantile(median, probs = .25, na.rm = TRUE),
            COI_q75 = quantile(median, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            LM_median = median(LM, na.rm=TRUE),
            LM_q25 = quantile(LM, probs = .25, na.rm = TRUE),
            LM_q75 = quantile(LM, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            DHA_median = median(DHA, na.rm=TRUE),
            DHA_q25 = quantile(DHA, probs = .25, na.rm = TRUE),
            DHA_q75 = quantile(DHA, probs = .75, na.rm = TRUE))


data %>% 
  summarize(N = n(),
            CQ_median = median(CQ, na.rm=TRUE),
            CQ_q25 = quantile(CQ, probs = .25, na.rm = TRUE),
            CQ_q75 = quantile(CQ, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            MDAQ_median = median(MDAQ, na.rm=TRUE),
            MDAQ_q25 = quantile(MDAQ, probs = .25, na.rm = TRUE),
            MDAQ_q75 = quantile(MDAQ, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            PQ_median = median(PQ, na.rm=TRUE),
            PQ_q25 = quantile(PQ, probs = .25, na.rm = TRUE),
            PQ_q75 = quantile(PQ, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            MQ_median = median(MQ, na.rm=TRUE),
            MQ_q25 = quantile(MQ, probs = .25, na.rm = TRUE),
            MQ_q75 = quantile(MQ, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            PND_median = median(PND, na.rm=TRUE),
            PND_q25 = quantile(PND, probs = .25, na.rm = TRUE),
            PND_q75 = quantile(PND, probs = .75, na.rm = TRUE))

data %>% 
  summarize(N = n(),
            PYR_median = median(PYR, na.rm=TRUE),
            PYR_q25 = quantile(PYR, probs = .25, na.rm = TRUE),
            PYR_q75 = quantile(PYR, probs = .75, na.rm = TRUE))

table(data$LM_cat)
table(data$CQ_cat)
table(data$MDAQ_cat)
table(data$PQ_cat)
table(data$MQ_cat)
table(data$PND_cat)
table(data$PYR_cat)



