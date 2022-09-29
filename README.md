# 2022_Tumwebaze_NatCom
Data sets and R code associated with Tumwebaze and Conrad et al (2022) Nature Comm

Data Files:
- uga_admbnda_ubos_20200824_shp.zip: Subset of Uganda administrative boundary shape files downloaded from Uganda Bureau of Statistics (adm2). Additional information about versions can be found in pdf contained in zip file. Used in generation of figure 1.
- PRX-04_k13_data.csv: PfK13 C469Y and A675V genotype calls for 2019 samples. Genotype calls made from MIP and dideoxy sequencing data available through NCBI (PRJNA655702 and MT857288-MT857721[https://www.ncbi.nlm.nih.gov/nuccore/?term=MT857288:MT857721[pacc]]. Used in generation of figure 1.
- exvivo_targetk13_MayAug21_wide.csv: contains data required to generate figures 2, 3, and 5, and Table 1, 2, and 3.  


Code Files:
- Fig1_Map.R: script used to generate "Figure 1: Map of Uganda". Requires extracted files in uga_admbnda_ubos_20200824_shp.zip and PRX-04_k13_data.csv
- Fig2_Fig3_Fig5.R: script used to generate figures 2, 3 and 5. Requires exvivo_targetk13_MayAug21_wide.csv
