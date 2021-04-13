#finding chemicals listed most often on FracFocus at sites near individuals in Fragile Families cohort

library(dplyr)
library(tidyr)

#file of FracFocus-reported chemicals linked to sites near individuals in Fragile Families cohort
chem <- read.csv("/home/HPA/HPA_group_datasets/Fracking/FracFocusRegistry_1-19cleaned.csv", stringsAsFactors=FALSE)

# filtering to ffid within 10 km (and chemicals linked to those people)
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")
dist_10 <- dist %>% filter(distance_ff_site <= 10)
ffid_10 <- unique(dist_10$ffid)

chem <- chem[chem$ffid %in% ffid_10,]

cas <- as.character(unique(chem$casnumber))

num <- c()

# getting counts for each cas number
for (i in 1:length(cas)){
cas_1 <- cas[i]
chem_filt <- chem %>% filter(casnumber == cas_1)
num_1 <- length(chem_filt$casnumber) # getting number of entries with this casnumber
num <- c(num, num_1)
}


df <- data.frame(cas, num)
