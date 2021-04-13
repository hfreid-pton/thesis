# analyzing SWPA EHP chemical data (find chemicals at higher concentration near fracking)
# file name on server: /home/HPA/Hannah/SWPA_chemanalysis_2_LS.R

library(dplyr)
library(tidyr)
library(ggplot2)

# reading chemical file
file <- read.csv("/home/HPA/Hannah/SWPA/SWPA_EHP_Summa_Data_2012-2019.csv", na.strings=c("", "","NA"),stringsAsFactors=F)
file[is.na(file)] = 0



# fixing unit issue
# there were only 4 entries in the ppbv unit that weren't just repeats
# of another entry that had the correct units (ug/m3)
# calculations done with online calculator
file$Units[file$Sample.ID == "18335-H763"] <- "ug/m3"
file$Units[file$Sample.ID == "18337-H760"] <- "ug/m3"
file$Units[file$Sample.ID == "18330-H759"] <- "ug/m3"
file$Units[file$Sample.ID == "18405-H761"] <- "ug/m3"
file$Units[file$Sample.ID == "18331-H762"] <- "ug/m3"
file$Acetone[file$Sample.ID == "18335-H763"] <- 2.9
file$Acetone[file$Sample.ID == "18337-H760"] <- 3.3
file$Acetone[file$Sample.ID == "18330-H759"] <- 2.6
file$Acetone[file$Sample.ID == "18405-H761"] <- 2.6
file$Acetone[file$Sample.ID == "18331-H762"] <- 2.6
file$Chloromethane[file$Sample.ID == "18335-H763"] <- 1.1
file$Chloromethane[file$Sample.ID == "18337-H760"] <- 1.1
file$Chloromethane[file$Sample.ID == "18330-H759"] <- 1.1
file$Chloromethane[file$Sample.ID == "18405-H761"] <- 1.2
file$Chloromethane[file$Sample.ID == "18331-H762"] <- 1.1
file$Dichlorodifluoromethane[file$Sample.ID == "18335-H763"] <- 2.7
file$Dichlorodifluoromethane[file$Sample.ID == "18337-H760"] <- 2.8
file$Dichlorodifluoromethane[file$Sample.ID == "18330-H759"] <- 2.8
file$Dichlorodifluoromethane[file$Sample.ID == "18405-H761"] <- 2.8
file$Dichlorodifluoromethane[file$Sample.ID == "18331-H762"] <- 2.8
file$Napthalene[file$Sample.ID == "18335-H763"] <- 3.7



#removing all entries with units of ppbv
file <- file %>% filter(Units == "ug/m3")

# remove cryogenic plants
file <- file %>% filter(Infrastucture != "cryogenic plant")


# making new baseline column
# 1=not baseline, 0=baseline
file$new_base <- ifelse(file$Baseline == "baseline", 0, 1)
file$new_base[is.na(file$new_base)] = 1


# seeing if chemical concentrations are signficantly different between baseline and non-baseline
# just comparing all non-baseline to all baseline

p_val <- c()
coef <- c()


for (i in 11:59){
#chemical <- chem[i]
model <- glm(file[,i] ~ file$new_base, data=file)
p_val1 <- coef(summary(model))[,4][2]
coef1 <- coef(summary(model))[,1][2]
p_val <- c(p_val, p_val1)
coef <- c(coef, coef1)
}

df_chem_analysis <- data.frame(chem, coef, p_val)

df_chem_exposure <- df_chem_analysis %>% filter(coef > 0)
df_chem_exposure <- df_chem_exposure %>% filter(p_val < 0.05)


# comparing non-baseline well pads to non-baseline compressor stations
file_type <- file %>% filter(new_base == 1)


# fixing infrastructure label
# compressor = 0, well-pad=1
file_type$new_infrastructure <- ifelse(file_type$Infrastucture == "compressor station", 0, 1)

p_val_type <- c()
coef_type <- c()

for (n in 11:58){
#chemical <- chem[i]
model2 <- glm(file_type[,n] ~ file_type$new_infrastructure, data=file_type)
p_val2 <- coef(summary(model2))[,4][2]
coef2 <- coef(summary(model2))[,1][2]
p_val_type <- c(p_val_type, p_val2)
coef_type <- c(coef_type, coef2)
}

df_chem_analysis_type <- data.frame(chem[1:48], coef_type, p_val_type)

df_chem_exposure_type <- df_chem_analysis_type %>% filter(p_val_type < 0.05)






















