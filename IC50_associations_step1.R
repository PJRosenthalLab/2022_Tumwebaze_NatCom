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

#setwd("") #recommend making directory specifically for IC50 analysis if RSA analysis is also planned
# Make new directory within working directory named "Wil_databases" within working directory where drug-specific databases will be written.

data_ic50_rsa <- read.csv("data_IC50_associations.csv")


data_LM <- subset(data_ic50_rsa, select = -c(MQ, CQ, PQ, DHA, PND, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_LM, "Wil_databases/MIP_LM_data.csv", row.names = FALSE)

data_MQ <- subset(data_ic50_rsa, select = -c(LM, CQ, PQ, DHA, PND, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_MQ, "Wil_databases/MIP_MQ_data.csv", row.names = FALSE)

data_CQ <- subset(data_ic50_rsa, select = -c(LM, MQ, PQ, DHA, PND, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_CQ, "Wil_databases/MIP_CQ_data.csv", row.names = FALSE)

data_PQ <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, DHA, PND, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_PQ, "Wil_databases/MIP_PQ_data.csv", row.names = FALSE)

data_DHA <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, PQ, PND, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_DHA, "Wil_databases/MIP_DHA_data.csv", row.names = FALSE)

data_PND <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, PQ, DHA, PYR, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_PND, "Wil_databases/MIP_PND_data.csv", row.names = FALSE)

data_PYR <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, PQ, DHA, PND, MDAQ, RSA_done, include, weighted_RSA))
write.csv(data_PYR, "Wil_databases/MIP_PYR_data.csv", row.names = FALSE)

data_MDAQ <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, PQ, DHA, PND, PYR, RSA_done, include, weighted_RSA))
write.csv(data_MDAQ, "Wil_databases/MIP_MDAQ_data.csv", row.names = FALSE)

data_RSA <- subset(data_ic50_rsa, select = -c(LM, MQ, CQ, PQ, DHA, PND, PYR, MDAQ, RSA_done))
write.csv(data_RSA, "Wil_databases/MIP_RSA_data.csv", row.names = FALSE)




