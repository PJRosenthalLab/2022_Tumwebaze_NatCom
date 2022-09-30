# 2022_Tumwebaze_NatCom
Data sets and R code associated with Tumwebaze and Conrad et al (2022) Nature Comm

Data Files:
- uga_admbnda_ubos_20200824_shp.zip: Subset of Uganda administrative boundary shape files downloaded from Uganda Bureau of Statistics (adm2). Additional information about versions can be found in pdf contained in zip file. Used in generation of figure 1.
- PRX-04_k13_data.csv: PfK13 C469Y and A675V genotype calls for 2019 samples. Genotype calls made from MIP and dideoxy sequencing data available through NCBI (PRJNA655702 and MT857288-MT857721[https://www.ncbi.nlm.nih.gov/nuccore/?term=MT857288:MT857721[pacc]]. Used in generation of figure 1.
- exvivo_targetk13_MayAug21_wide.csv: contains data required to generate figures 2, 3, and 5, and Table 1 (Isolates with RSA and/or IC50 results), 2, and 3. 
- exvivo_all_wide.csv: contains data required to populate table 1 (Isolates with only IC50 results)
- artR_genes_wide.csv: contains data required to generate Figure 4 and perform genotype-phenotype associations for RSA data.
- data_IC50_associations.csv: contains data required to populate Table 4 and perform genotype-phenotype associations for IC50 data.
- Supplemental_Table_9_Design of drug resistance MIP panel.xlsx: contains details on probe design for drug resistance MIP panel

Code Files:
- Fig1_Map.R: script used to generate "Figure 1: Map of Uganda". Requires extracted files in uga_admbnda_ubos_20200824_shp.zip and PRX-04_k13_data.csv
- Fig2_Fig3_Fig5.R: script used to generate figures 2, 3 and 5. Requires exvivo_targetk13_MayAug21_wide.csv
- Tables1-3.R: script used to populate Tables1-3 (for table 1 limited do Isolates with RSA and/or IC50 results). Requires exvivo_targetk13_MayAug21_wide.csv 
- Tables1_all_IC50: script used to poulate Table 1 (Isolates with only IC50 results). Requires exvivo_all_wide.csv
- RSA_associations.R and RSA_associations_step2.R: used for genotype-phenotype associations for RSA data and to generate dataframe to make Figure 4. Requires artR_genes_wide.csv data.
- IC50_associations_step1.R and RSA_associations_step2.R: used for genotype-phenotype associations for IC50 data and to generate dataframe to make table 4. Requires data_IC50_associations.csv
