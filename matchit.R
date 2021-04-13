#matching controls to exposed individuals
#server code name: /home/HPA/Hannah/matchit_11-1.R 

library(dplyr)
library(MatchIt)
library(tidyr)
library(stringr)

# reading file of distances to fracking sites
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")

# defining who is exposed
dist_10 <- dist %>% filter(distance_ff_site <= 10)
all_ffid <- unique(dist$ffid)

ffid_10 <- dist_10$ffid
counts_10 <- c()

for (i in 1:length(all_ffid)){
num_10 <- length(ffid_10[ffid_10 == all_ffid[i]])
counts_10 <- append(counts_10, num_10)
}

num_wells <- data.frame(all_ffid, counts_10)

# changing from ffid to ffFam ID
ffam_test <- str_sub(num_wells$all_ffid, end=-5)
num_wells$ffFam <- ffam_test

# reading in survey data--need matching variables
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

exposed <- num_wells %>% filter(counts_10 > 0)


#finding potential controls to be used for matching
within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam


num_wells_10 <- num_wells
num_wells_10$exposure <- ifelse(num_wells_10$counts_10 > 0, 1, 0)
num_wells_10 <- num_wells_10[(num_wells_10$ffFam %in% controls_ffFam) | (num_wells_10$ffFam %in% exposed$ffFam),]

#keeping only exposed and potential controls in survey data
merged2 <- surv[surv$ffFam %in% num_wells_10$ffFam,]

merged3 <- merge(merged2, num_wells_10, by="ffFam")


#dropping NAs/error codes for variables used in matching
merged4 <- drop_na(merged3,c("cm1ethrace", "cm1bsex", "cm5povca", "p5h15c", "m1g4", "k6d46", "ck6pcgrel", "cp6povca"))


merged4 <- merged4[merged4$cm1ethrace > 0,]
merged4 <- merged4[merged4$cm1ethrace < 4,]
merged4 <- merged4[merged4$cm1bsex > 0,]
merged4 <- merged4[merged4$cm5povca > 0,]
merged4 <- merged4[merged4$p5h15c >= 0,]
merged4 <- merged4[merged4$m1g4 > 0,]
merged4 <- merged4[merged4$k6d46 >= 0,]
merged4 <- merged4[merged4$ck6pcgrel > 0,]
merged4 <- merged4[merged4$cp6povca > 0,]

#keeping only people who didn't move or moved but didn't switch schools
merged4 <- merged4 %>% filter(p6j1 == 2 | p6j5 == 0)


merged4 <- merged4 %>% select(ffFam, cm1ethrace,exposure,cm1bsex, cm5povca, p5h15c, m1g4,k6d46, ck6pcgrel, cp6povca)





#making cm1ethrace into dummary variable(white=reference)
merged4$his <- ifelse(merged4$cm1ethrace == 3, 1, 0)
merged4$afr <- ifelse(merged4$cm1ethrace == 2, 1, 0)

# matching
m.out <- matchit(exposure ~ his + afr + cm1bsex + cm5povca + p5h15c + m1g4 + k6d46 + ck6pcgrel + cp6povca, data=merged4, method="nearest", ratio=1)







