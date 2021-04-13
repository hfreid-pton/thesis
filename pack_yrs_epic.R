#looking for association between PACKYRs score and self-reported smoking measures

library(dplyr)
library(tidyr)
library(stringr)

#new epigenetic clock data (not bmiq normalized)
clocks <- read.csv("/home/HPA/HPA_group_datasets/MethylationData/January2021/epic_epigeneticAges.csv")

# removing samples with potential issues
prob1 <- readRDS("/home/HPA/Hannah/BenCode/datafiles/probIDsEpic.RDS")
prob2 <- read.csv("/home/HPA/Hannah/prob_allEPIC_hannah.csv", stringsAsFactors=FALSE)
prob3 <- read.csv("/home/HPA/Hannah/prob_clock_allepic_hannah.csv", stringsAsFactors=FALSE)
prob <- c(prob1, prob2$prob_both, prob3$x)
prob <- unique(prob)


# loading survey data which has self-reported smoking measures
surv2 <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")
surv2 <- surv2 %>% select(ffFam, cm1bsex)
colnames(surv2) <- c("ffFam", "cm1bsex2")

# separating ffFam from wave in SampleID
ffFam <- str_sub(clocks$ffID, end=-5)
clocks$ffFam <- ffFam

#removing problem samples
clocks <- clocks[!clocks$ffFam %in% prob,]

clocks <- merge(clocks, surv2, by="ffFam")

# filtering gender discrepancies
clocks_male <- clocks %>% filter(predictedGender == "male" & cm1bsex2==1)
clocks_female <- clocks %>% filter(predictedGender == "female" & cm1bsex2 ==2)
clocks <- rbind(clocks_male, clocks_female)


clocks_9 <- clocks %>% filter(Age < 13)
clocks_15 <- clocks %>% filter(Age > 13)
ffFam_both <- intersect(clocks_15$ffFam, clocks_9$ffFam)
clocks_both <- clocks[clocks$ffFam %in% ffFam_both,]



#PackYRS age adj residual
PACK_res <- residuals(lm(DNAmPACKYRS ~ Age, data=clocks_both), na.action=na.exclude)
clocks_both$PACK_res <- PACK_res


# making clock df with age9 and age15 scores in same row for same person (i.e. one row/person) 
unique_ffFam_both <- unique(clocks_both$ffFam)
age_9_both <- c()
age_15_both <- c()
age_9_pack_res <- c()
age_15_pack_res <- c()


clocks_9_new <- clocks_both %>% filter(Age < 13)
clocks_15_new <- clocks_both %>% filter(Age > 13)
for (i in 1:length(unique_ffFam_both)){
filt_9 <- clocks_9_new %>% filter(ffFam == unique_ffFam_both[i])
filt_15 <- clocks_15_new %>% filter(ffFam == unique_ffFam_both[i])
age_9_pack_res <- c(age_9_pack_res, filt_9$PACK_res[1])
age_15_pack_res <- c(age_15_pack_res, filt_15$PACK_res[1])
}


ffFam <- unique_ffFam_both
df_both <- data.frame(ffFam, age_9_pack_res, age_15_pack_res)


# merging with survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")
surv <- drop_na(surv, c("m1g4", "k6d46", "p5h15c"))
surv <- surv[surv$m1g4 > 0,]
surv <- surv[surv$p5h15c >= 0,]
surv <- surv[surv$k6d46 >= 0,]

#fixing smoking data coding
surv$m1g4_corrected1 <- ifelse(surv$m1g4 == 4, 0, surv$m1g4)
surv$m1g4_corrected2 <- ifelse(surv$m1g4_corrected1 == 3, 4, surv$m1g4_corrected1)
surv$m1g4_corrected3 <- ifelse(surv$m1g4_corrected2 == 1, 3, surv$m1g4_corrected2)
surv$m1g4_corrected <- ifelse(surv$m1g4_corrected3 == 4, 1, surv$m1g4_corrected3)
surv$k6d46_corrected <- ifelse(surv$k6d46 == 1, 0, surv$k6d46)

df_both_merged <- merge(df_both, surv)
write.csv(df_both_merged, "/home/HPA/Hannah/PACKYRs_mergedsmoking_2-20.csv")

# looking for associations
model_15_base <- lm(age_15_pack_res ~ m1g4_corrected, data=df_both_merged)
model_15_base_9 <- lm(age_15_pack_res ~ m1g4_corrected + p5h15c, data=df_both_merged)
model_15_base_9_15 <- lm(age_15_pack_res ~ m1g4_corrected + p5h15c + k6d46_corrected, data=df_both_merged)

model_9_base <- lm(age_9_pack_res ~ m1g4_corrected, data=df_both_merged)
model_9_base_9 <- lm(age_9_pack_res ~ m1g4_corrected + p5h15c, data=df_both_merged)




