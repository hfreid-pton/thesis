#descriptive stats

library(dplyr)
library(tidyr)
library(stringr)

# loading distance to fracking site data and finding who is exposed and who is potential control
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")
dist_10 <- dist %>% filter(distance_ff_site <= 10)
all_ffid <- unique(dist$ffid)


ffid_10 <- dist_10$ffid
counts_10 <- c()

for (i in 1:length(all_ffid)){
num_10 <- length(ffid_10[ffid_10 == all_ffid[i]])
counts_10 <- append(counts_10, num_10)
}

num_wells <- data.frame(all_ffid, counts_10)

# changing from ffid to ffFam id
ffam_test <- str_sub(num_wells$all_ffid, end=-5)
num_wells$ffFam <- ffam_test

# loading survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

#defining exposed vs. potential controls
exposed <- num_wells %>% filter(counts_10 > 0)


within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam


num_wells_10 <- num_wells
num_wells_10$exposure <- ifelse(num_wells_10$counts_10 > 0, 1, 0)
num_wells_10 <- num_wells_10[(num_wells_10$ffFam %in% controls_ffFam) | (num_wells_10$ffFam %in% exposed$ffFam),]

#keeping just survey data for people in exposed/potential control groups
surv2 <- surv[surv$ffFam %in% num_wells_10$ffFam,]
merged <- merge(num_wells_10, surv2, by="ffFam")

#getting rid of NAs/error codes
merged2 <- drop_na(merged,c("cm1ethrace", "cm1bsex", "cm5povca", "p5h15c", "m1g4", "k6d46", "ck6pcgrel", "cp6povca"))

merged2 <- merged2[merged2$cm1ethrace > 0,]
merged2 <- merged2[merged2$cm1ethrace < 4,]
merged2 <- merged2[merged2$cm1bsex > 0,]
merged2 <- merged2[merged2$cm5povca > 0,]
merged2 <- merged2[merged2$p5h15c >= 0,]
merged2 <- merged2[merged2$m1g4 > 0,]
merged2 <- merged2[merged2$k6d46 >= 0,]
merged2 <- merged2[merged2$ck6pcgrel > 0,]
merged2 <- merged2[merged2$cp6povca > 0,]

#keeping only people who did not move or who moved but did not switch schools
merged2 <- merged2 %>% filter(p6j1 == 2 | p6j5 == 0)

#changing poverty to binary
merged2$pov5 <- ifelse(merged2$cm5povca > 2, 0, 1)
merged2$pov6 <- ifelse(merged2$cp6povca > 2, 0, 1)

#changing smoking exposure to binary
merged2$smoke1 <- ifelse(merged2$m1g4 < 4, 1, 0)
merged2$smoke2 <- ifelse(merged2$p5h15c > 0, 1, 0)
merged2$smoke3 <- ifelse(merged2$k6d46 > 1, 1, 0)


merged2_exp <- merged2 %>% filter(exposure == 1)
merged2_c <- merged2 %>% filter(exposure == 0)

# finding number of indivduals in each group for each variable
# did calculation to go from number to percent offline on calculator by hand
table(merged2_exp$cm1ethrace)
table(merged2_c$cm1ethrace)
table(merged2_exp$cm1bsex)
table(merged2_c$cm1bsex)
table(merged2_exp$pov5)
table(merged2_c$pov5)
table(merged2_exp$pov6)
table(merged2_c$pov6)
table(merged2_exp$ck6pcgrel)
table(merged2_c$ck6pcgrel)
table(merged2_exp$smoke1)
table(merged2_c$smoke1)
table(merged2_exp$smoke2)
table(merged2_c$smoke2)
table(merged2_exp$smoke3)
table(merged2_c$smoke3)










