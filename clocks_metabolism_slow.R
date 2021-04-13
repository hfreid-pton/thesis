# comparing methylation age clocks in exposed slow metabolizers with control slow metabolizers
# server code name: /home/HPA/Hannah/clocks_metabolism_slow.R

library(dplyr)
library(tidyr)
library(stringr)


# read in metabolism scores and create binary scores
score <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/9var_scores_3-14.csv")
score$ffFam <- str_sub(score$Indiv, end=-5)
score$binary <- ifelse(score$total < 5,0,1)
score$total_alc <- (score$score_2e16 + score$score_adh1c + score$score_adh1b + score$score_e15)
score$binary_alc <- ifelse(score$total_alc < 2,0,1)



# read in methylation clock data
clocks <- read.csv("/home/HPA/HPA_group_datasets/MethylationData/January2021/epic_epigeneticAges.csv")

# read in list of problem samples
prob1 <- readRDS("/home/HPA/Hannah/BenCode/datafiles/probIDsEpic.RDS")
prob2 <- read.csv("/home/HPA/Hannah/prob_allEPIC_hannah.csv", stringsAsFactors=FALSE)
prob3 <- read.csv("/home/HPA/Hannah/prob_clock_allepic_hannah.csv", stringsAsFactors=FALSE)
prob <- c(prob1, prob2$prob_both, prob3$x)
prob <- unique(prob)


# separating ffFam from wave in SampleID
ffFam <- str_sub(clocks$ffID, end=-5)
clocks$ffFam <- ffFam


#removing problem samples
clocks <- clocks[!clocks$ffFam %in% prob,]


# keeping only people with scores and clock data
clocks <- clocks[clocks$ffFam %in% score$ffFam,]


#keeping only people with slow metabolism scores
score_slow <- score$ffFam[score$binary == 1] #this is changed to score$binary_alc == 1 for using alc score
clocks <- clocks[clocks$ffFam %in% score_slow,]


# read in survey data for filtering gender discrepancies
surv2 <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")
surv2 <- surv2 %>% select(ffFam, cm1bsex)
colnames(surv2) <- c("ffFam", "cm1bsex2")
clocks <- merge(clocks, surv2, by="ffFam")

# filtering gender discrepancies
clocks_male <- clocks %>% filter(predictedGender == "male" & cm1bsex2==1)
clocks_female <- clocks %>% filter(predictedGender == "female" & cm1bsex2 ==2)
clocks <- rbind(clocks_male, clocks_female)


# keeping only people with data at both ages
clocks_9 <- clocks %>% filter(Age < 13)
clocks_15 <- clocks %>% filter(Age > 13)
clocks_9 <- clocks_9[complete.cases(clocks_9),]
clocks_15 <- clocks_15[complete.cases(clocks_15),]
ffFam_both <- intersect(clocks_15$ffFam, clocks_9$ffFam)
clocks_both <- clocks[clocks$ffFam %in% ffFam_both,]


#calculating age adj residual for each clock (some clocks already
# have this calculated--by Lisa Schneper)

#PackYRS
PACK_res <- residuals(lm(DNAmPACKYRS ~ Age, data=clocks_both, na.action=na.exclude))
clocks_both$PACK_res <- PACK_res


#DunedinPoAm_38
Dunedin_res <- residuals(lm(DunedinPoAm_38 ~ Age, data=clocks_both, na.action=na.exclude))
clocks_both$Dunedin_res <- Dunedin_res

#Hannum Age
Hannum_res <- residuals(lm(DNAmHannumAge ~ Age, data=clocks_both, na.action=na.exclude))
clocks_both$Hannum_res <- Hannum_res

#Ped Age
Ped_res <- residuals(lm(PedBE_age ~ Age, data=clocks_both, na.action=na.exclude))
clocks_both$Ped_res <- Ped_res

#remaking df to have age 9 and age 15 meth ages in same row (i.e. one row per person)
unique_ffFam_both <- unique(clocks_both$ffFam)
age_9_both <- c()
age_15_both <- c()
chr_age_9 <- c()
chr_age_15 <- c()
age_9_pack <- c()
age_15_pack <- c()
age_9_pack_res <- c()
age_15_pack_res <- c()
age_9_pheno <- c()
age_15_pheno <- c()
age_9_pheno_res <- c()
age_15_pheno_res <- c()
age_9_grim <- c()
age_15_grim <- c()
age_9_grim_res <- c()
age_15_grim_res <- c()
age_9_dunedin <- c()
age_15_dunedin <- c()
age_9_dunedin_res <- c()
age_15_dunedin_res <- c()
age_9_dnam <- c()
age_15_dnam <- c()
age_9_dnam_res <- c()
age_15_dnam_res <- c()
age_9_hannum <- c()
age_15_hannum <- c()
age_9_hannum_res <- c()
age_15_hannum_res <- c()
age_9_ped <- c()
age_15_ped <- c()
age_9_ped_res <- c()
age_15_ped_res <- c()
age_9_bioage4 <- c()
age_15_bioage4 <- c()
age_9_bioage4_res <- c()
age_15_bioage4_res <- c()

clocks_9_new <- clocks_both %>% filter(Age < 13)
clocks_15_new <- clocks_both %>% filter(Age > 13)

for (i in 1:length(unique_ffFam_both)){
filt_9 <- clocks_9_new %>% filter(ffFam == unique_ffFam_both[i])
filt_15 <- clocks_15_new %>% filter(ffFam == unique_ffFam_both[i])
chr_age_9 <- c(chr_age_9, filt_9$Age[1])
chr_age_15 <- c(chr_age_15, filt_15$Age[1])
age_9_pack <- c(age_9_pack, filt_9$DNAmPACKYRS[1])
age_15_pack <- c(age_15_pack, filt_15$DNAmPACKYRS[1])
age_9_pack_res <- c(age_9_pack_res, filt_9$PACK_res[1])
age_15_pack_res <- c(age_15_pack_res, filt_15$PACK_res[1]) 
age_9_pheno <- c(age_9_pheno, filt_9$DNAm_phenoAge[1])
age_15_pheno <- c(age_15_pheno, filt_15$DNAm_phenoAge[1])
age_9_pheno_res <- c(age_9_pheno_res, filt_9$AgeAccelPheno[1])
age_15_pheno_res <- c(age_15_pheno_res, filt_15$AgeAccelPheno[1])
age_9_grim <- c(age_9_grim, filt_9$DNAmGrimAge[1])
age_15_grim <- c(age_15_grim, filt_15$DNAmGrimAge[1])
age_9_grim_res <- c(age_9_grim_res, filt_9$AgeAccelGrim[1])
age_15_grim_res <- c(age_15_grim_res, filt_15$AgeAccelGrim[1])
age_9_dunedin <- c(age_9_dunedin, filt_9$DunedinPoAm_38[1])
age_15_dunedin <- c(age_15_dunedin, filt_15$DunedinPoAm_38[1])
age_9_dunedin_res <- c(age_9_dunedin_res, filt_9$Dunedin_res[1])
age_15_dunedin_res <- c(age_15_dunedin_res, filt_15$Dunedin_res[1])
age_9_dnam <- c(age_9_dnam, filt_9$DNAmAge[1])
age_15_dnam <- c(age_15_dnam, filt_15$DNAmAge[1])
age_9_dnam_res <- c(age_9_dnam_res, filt_9$IEAA[1])
age_15_dnam_res <- c(age_15_dnam_res, filt_15$IEAA[1])
age_9_ped <- c(age_9_ped, filt_9$PedBE_age[1])
age_15_ped <- c(age_15_ped, filt_15$PedBE_age[1])
age_9_ped_res <- c(age_9_ped_res, filt_9$Ped_res[1])
age_15_ped_res <- c(age_15_ped_res, filt_15$Ped_res[1])
age_9_hannum <- c(age_9_hannum, filt_9$DNAmHannumAge[1])
age_15_hannum <- c(age_15_hannum, filt_15$DNAmHannumAge[1])
age_9_hannum_res <- c(age_9_hannum_res, filt_9$Hannum_res[1])
age_15_hannum_res <- c(age_15_hannum_res, filt_15$Hannum_res[1])
age_9_dnampai <- c(age_9_dnampai, filt_9$DNAmpai_1[1])
age_15_dnampai <- c(age_15_dnampai, filt_15$DNAmpai_1[1])
age_9_dnampai_res <- c(age_9_dnampai_res, filt_9$Adj.DNAmPAI_1[1])
age_15_dnampai_res <- c(age_15_dnampai_res, filt_15$Adj.DNAmPAI_1[1])
age_9_bioage4 <- c(age_9_bioage4, filt_9$BioAge4HAStatic[1])
age_15_bioage4 <- c(age_15_bioage4, filt_15$BioAge4HAStatic[1])
age_9_bioage4_res <- c(age_9_bioage4_res, filt_9$EEAA[1])
age_15_bioage4_res <- c(age_15_bioage4_res, filt_15$EEAA[1])
}

ffFam <- unique_ffFam_both

df_both <- data.frame(ffFam, chr_age_9, chr_age_15, age_9_pack, age_15_pack, age_9_pack_res, age_15_pack_res, age_9_pheno, age_15_pheno, age_9_pheno_res, age_15_pheno_res, age_9_grim, age_15_grim, age_9_grim_res, age_15_grim_res, age_9_dunedin, age_15_dunedin, age_9_dunedin_res, age_15_dunedin_res, age_9_dnam, age_15_dnam, age_9_dnam_res, age_15_dnam_res, age_9_ped, age_15_ped, age_9_ped_res, age_15_ped_res, age_9_hannum, age_15_hannum, age_9_hannum_res, age_15_hannum_res, age_9_dnampai, age_15_dnampai, age_9_dnampai_res, age_15_dnampai_res, age_9_bioage4, age_15_bioage4, age_9_bioage4_res, age_15_bioage4_res)

df_both$chrage_diff <- chr_age_15 - chr_age_9
chrage_diff <- df_both$chrage_diff

# adding columns for change in age adj residual and pace of epigenetic aging
df_both <- df_both %>% mutate(pack_dif=(age_15_pack_res - age_9_pack_res))
df_both <- df_both %>% mutate(pack_pace = ((age_15_pack - age_9_pack)/chrage_diff))
df_both <- df_both %>% mutate(pheno_dif = (age_15_pheno_res - age_9_pheno_res))
df_both <- df_both %>% mutate(pheno_pace = ((age_15_pheno-age_9_pheno)/chrage_diff))
df_both <- df_both %>% mutate(grim_dif = (age_15_grim_res - age_9_pheno_res))
df_both <- df_both %>% mutate(grim_pace = ((age_15_grim - age_9_grim)/chrage_diff))
df_both <- df_both %>% mutate(dunedin_dif = (age_15_dunedin_res - age_9_dunedin_res))
df_both <- df_both %>% mutate(dunedin_pace=((age_15_dunedin - age_9_dunedin)/chrage_diff))
df_both <- df_both %>% mutate(dnam_dif = (age_15_dnam_res - age_9_dnam_res))
df_both <- df_both %>% mutate(dnam_pace= ((age_15_dnam - age_9_dnam)/chrage_diff))
df_both <- df_both %>% mutate(ped_dif = (age_15_ped_res - age_9_ped_res))
df_both <- df_both %>% mutate(ped_pace = ((age_15_ped - age_9_ped)/chrage_diff))
df_both <- df_both %>% mutate(hannum_dif = (age_15_hannum_res - age_9_hannum_res))
df_both <- df_both %>% mutate(hannum_pace = ((age_15_hannum - age_9_hannum)/chrage_diff))
df_both <- df_both %>% mutate(bioage4_dif = (age_15_bioage4_res - age_9_bioage4_res))
df_both <- df_both %>% mutate(bioage4_pace = ((age_15_bioage4 - age_9_bioage4)/chrage_diff))



# setting up vectors to hold p-values, coefficients, and means before running 100 iterations
p_pack_15res <- c()
coef_pack_15res <- c()
p_pack_9res <- c()
coef_pack_9res <- c()
p_pack_res <- c()
coef_pack_res <- c()
p_pack_pace <- c()
coef_pack_pace <- c()
mean_pack_15res_exp <- c()
mean_pack_15res_c <- c()
mean_pack_res_exp <- c()
mean_pack_res_c <- c()
mean_pack_pace_exp <- c()
mean_pack_pace_c <- c()


p_pheno_15res <- c()
coef_pheno_15res <- c()
p_pheno_9res <- c()
coef_pheno_9res <- c()
p_pheno_res <- c()
coef_pheno_res <- c()
p_pheno_pace <- c()
coef_pheno_pace <-c()
mean_pheno_15res_exp <- c()
mean_pheno_15res_c <- c()
mean_pheno_res_exp <- c()
mean_pheno_res_c <- c()
mean_pheno_pace_exp <- c()
mean_pheno_pace_c <- c()



p_grim_15res <- c()
coef_grim_15res <- c()
p_grim_9res <- c()
coef_grim_9res <- c()
p_grim_res <- c()
coef_grim_res <- c()
p_grim_pace <- c()
coef_grim_pace <- c()
mean_grim_15res_exp <- c()
mean_grim_15res_c <- c()
mean_grim_res_exp <- c()
mean_grim_res_c <- c()
mean_grim_pace_exp <- c()
mean_grim_pace_c <- c()




p_dunedin_15res <- c()
coef_dunedin_15res <- c()
p_dunedin_9res <- c()
coef_dunedin_9res <- c()
p_dunedin_res <- c()
coef_dunedin_res <- c()
p_dunedin_pace <- c()
coef_dunedin_pace <- c()
mean_dunedin_15res_exp <- c()
mean_dunedin_15res_c <- c()
mean_dunedin_res_exp <- c()
mean_dunedin_res_c <- c()
mean_dunedin_pace_exp <- c()
mean_dunedin_pace_c <- c()



p_dnam_15res <- c()
coef_dnam_15res <- c()
p_dnam_9res <- c()
coef_dnam_9res <- c()
p_dnam_res <- c()
coef_dnam_res <- c()
p_dnam_pace <- c()
coef_dnam_pace <- c()
mean_dnam_15res_exp <- c()
mean_dnam_15res_c <- c()
mean_dnam_res_exp <- c()
mean_dnam_res_c <- c()
mean_dnam_pace_exp <- c()
mean_dnam_pace_c <- c()



p_ped_15res <- c()
coef_ped_15res <- c()
p_ped_9res <- c()
coef_ped_9res <- c()
p_ped_res <- c()
coef_ped_res <- c()
p_ped_pace <- c()
coef_ped_pace <- c()
mean_ped_15res_exp <- c()
mean_ped_15res_c <- c()
mean_ped_res_exp <- c()
mean_ped_res_c <- c()
mean_ped_pace_exp <- c()
mean_ped_pace_c <- c()



p_hannum_15res <- c()
coef_hannum_15res <- c()
p_hannum_9res <- c()
coef_hannum_9res <- c()
p_hannum_res <- c()
coef_hannum_res <- c()
p_hannum_pace <- c()
coef_hannum_pace <-c()
mean_hannum_15res_exp <- c()
mean_hannum_15res_c <- c()
mean_hannum_res_exp <- c()
mean_hannum_res_c <- c()
mean_hannum_pace_exp <- c()
mean_hannum_pace_c <- c()

p_bioage4_15res <- c()
coef_bioage4_15res <- c()
p_bioage4_9res <- c()
coef_bioage4_9res <- c()
p_bioage4_res <- c()
coef_bioage4_res <- c()
p_bioage4_pace <- c()
coef_bioage4_pace <- c()
mean_bioage4_15res_exp <- c()
mean_bioage4_15res_c <- c()
mean_bioage4_res_exp <- c()
mean_bioage4_res_c <- c()
mean_bioage4_pace_exp <- c()
mean_bioage4_pace_c <- c()



# loop to compare exposed with 100 iterations of matched controls

for (r in 1:100){
print(r)

#matching
#metabolism specific matching code just makes sure only individuals with metaboslism scores 
# (and scores that fit criteria if slow only comaprison) are included for matching

source("/home/HPA/Hannah/matchit_11-1_metabolism.R")
match_df <- match.data(m.out)
merged_both_nn <- merge(df_both, match_df, by="ffFam")

# seeing if exposure is associated with methylation age (for each clock)
model_nn_15_pack <- glm(age_15_pack_res ~ exposure, data=merged_both_nn)
model_nn_9_pack <- glm(age_9_pack_res ~ exposure, data=merged_both_nn)
model_nn_dif_pack <- glm(pack_dif ~ exposure, data=merged_both_nn)
model_nn_pace_pack <- glm(pack_pace ~ exposure, data=merged_both_nn)

model_nn_15_pheno <- glm(age_15_pheno_res ~ exposure, data=merged_both_nn)
model_nn_9_pheno <- glm(age_9_pheno_res ~ exposure, data=merged_both_nn)
model_nn_dif_pheno <- glm(pheno_dif ~ exposure, data=merged_both_nn)
model_nn_pace_pheno <- glm(pheno_pace ~ exposure, data=merged_both_nn)

model_nn_15_grim <- glm(age_15_grim_res ~ exposure, data=merged_both_nn)
model_nn_9_grim <- glm(age_9_grim_res ~ exposure, data=merged_both_nn)
model_nn_dif_grim <- glm(grim_dif ~ exposure, data=merged_both_nn)
model_nn_pace_grim <- glm(grim_pace ~ exposure, data=merged_both_nn)

model_nn_15_dunedin <- glm(age_15_dunedin_res ~ exposure, data=merged_both_nn)
model_nn_9_dunedin <- glm(age_9_dunedin_res ~ exposure, data=merged_both_nn)
model_nn_dif_dunedin <- glm(dunedin_dif ~ exposure, data=merged_both_nn)
model_nn_pace_dunedin <- glm(dunedin_pace ~ exposure, data=merged_both_nn)

model_nn_15_dnam <- glm(age_15_dnam_res ~ exposure, data=merged_both_nn)
model_nn_9_dnam <- glm(age_9_dnam_res ~ exposure, data=merged_both_nn)
model_nn_dif_dnam <- glm(dnam_dif ~ exposure, data=merged_both_nn)
model_nn_pace_dnam <- glm(dnam_pace ~ exposure, data=merged_both_nn)

model_nn_15_ped <- glm(age_15_ped_res ~ exposure, data=merged_both_nn)
model_nn_9_ped <- glm(age_9_ped_res ~ exposure, data=merged_both_nn)
model_nn_dif_ped <- glm(ped_dif ~ exposure, data=merged_both_nn)
model_nn_pace_ped <- glm(ped_pace ~ exposure, data=merged_both_nn)

model_nn_15_hannum <- glm(age_15_hannum_res ~ exposure, data=merged_both_nn)
model_nn_9_hannum <- glm(age_9_hannum_res ~ exposure, data=merged_both_nn)
model_nn_dif_hannum <- glm(hannum_dif ~ exposure, data=merged_both_nn)
model_nn_pace_hannum <- glm(hannum_pace ~ exposure, data=merged_both_nn)

model_nn_15_bioage4 <- glm(age_15_bioage4_res ~ exposure, data=merged_both_nn)
model_nn_9_bioage4 <- glm(age_9_bioage4_res ~ exposure, data=merged_both_nn)
model_nn_dif_bioage4 <- glm(bioage4_dif ~ exposure, data=merged_both_nn)
model_nn_pace_bioage4 <- glm(bioage4_pace ~ exposure, data=merged_both_nn)

# dfs to be used for finding mean methylation ages for each group for each clock
merged_both_nn_exp <- merged_both_nn %>% filter(exposure == 1)
merged_both_nn_c <- merged_both_nn %>% filter(exposure == 0)


# saving coefficients/pvals from each iteration of the loop
p_pack_15res_1 <- coef(summary(model_nn_15_pack))[,4][2]
p_pack_15res <- c(p_pack_15res, p_pack_15res_1)
coef_pack_15res_1 <- coef(summary(model_nn_15_pack))[,1][2]
coef_pack_15res <- c(coef_pack_15res, coef_pack_15res_1)

p_pack_9res_1 <- coef(summary(model_nn_9_pack))[,4][2]
p_pack_9res <- c(p_pack_9res, p_pack_9res_1)
coef_pack_9res_1 <- coef(summary(model_nn_9_pack))[,1][2]
coef_pack_9res <- c(coef_pack_9res, coef_pack_9res_1)

p_pack_res_1 <- coef(summary(model_nn_dif_pack))[,4][2]
p_pack_res <- c(p_pack_res, p_pack_res_1)
coef_pack_res_1 <- coef(summary(model_nn_dif_pack))[,1][2]
coef_pack_res <- c(coef_pack_res, coef_pack_res_1)

p_pack_pace_1 <- coef(summary(model_nn_pace_pack))[,4][2]
p_pack_pace <- c(p_pack_pace, p_pack_pace_1)
coef_pack_pace_1 <- coef(summary(model_nn_pace_pack))[,1][2]
coef_pack_pace <- c(coef_pack_pace, coef_pack_pace_1)

#saving mean for each clock for each group from each iteration 
# (this and above to save coeff/pvals repeats for each clock below)
mean_pack_15res_exp <- c(mean_pack_15res_exp, mean(merged_both_nn_exp$age_15_pack_res))
mean_pack_15res_c <- c(mean_pack_15res_c, mean(merged_both_nn_c$age_15_pack_res))
mean_pack_res_exp <- c(mean_pack_res_exp, mean(merged_both_nn_exp$pack_dif))
mean_pack_res_c <- c(mean_pack_res_c, mean(merged_both_nn_c$pack_dif))
mean_pack_pace_exp <- c(mean_pack_pace_exp, mean(merged_both_nn_exp$pack_pace))
mean_pack_pace_c <- c(mean_pack_pace_c, mean(merged_both_nn_c$pack_pace))


p_pheno_15res_1 <- coef(summary(model_nn_15_pheno))[,4][2]
p_pheno_15res <- c(p_pheno_15res, p_pheno_15res_1)
coef_pheno_15res_1 <- coef(summary(model_nn_15_pheno))[,1][2]
coef_pheno_15res <- c(coef_pheno_15res, coef_pheno_15res_1)

p_pheno_9res_1 <- coef(summary(model_nn_9_pheno))[,4][2]
p_pheno_9res <- c(p_pheno_9res, p_pheno_9res_1)
coef_pheno_9res_1 <- coef(summary(model_nn_9_pheno))[,1][2]
coef_pheno_9res <- c(coef_pheno_9res, coef_pheno_9res_1)

p_pheno_res_1 <- coef(summary(model_nn_dif_pheno))[,4][2]
p_pheno_res <- c(p_pheno_res, p_pheno_res_1)
coef_pheno_res_1 <- coef(summary(model_nn_dif_pheno))[,1][2]
coef_pheno_res <- c(coef_pheno_res, coef_pheno_res_1)

p_pheno_pace_1 <- coef(summary(model_nn_pace_pheno))[,4][2]
p_pheno_pace <- c(p_pheno_pace, p_pheno_pace_1)
coef_pheno_pace_1 <- coef(summary(model_nn_pace_pheno))[,1][2]
coef_pheno_pace <- c(coef_pheno_pace, coef_pheno_pace_1)

mean_pheno_15res_exp <- c(mean_pheno_15res_exp, mean(merged_both_nn_exp$age_15_pheno_res))
mean_pheno_15res_c <- c(mean_pheno_15res_c, mean(merged_both_nn_c$age_15_pheno_res))
mean_pheno_res_exp <- c(mean_pheno_res_exp, mean(merged_both_nn_exp$pheno_dif))
mean_pheno_res_c <- c(mean_pheno_res_c, mean(merged_both_nn_c$pheno_dif))
mean_pheno_pace_exp <- c(mean_pheno_pace_exp, mean(merged_both_nn_exp$pheno_pace))
mean_pheno_pace_c <- c(mean_pheno_pace_c, mean(merged_both_nn_c$pheno_pace))




p_grim_15res_1 <- coef(summary(model_nn_15_grim))[,4][2]
p_grim_15res <- c(p_grim_15res, p_grim_15res_1)
coef_grim_15res_1 <- coef(summary(model_nn_15_grim))[,1][2]
coef_grim_15res <- c(coef_grim_15res, coef_grim_15res_1)

p_grim_9res_1 <- coef(summary(model_nn_9_grim))[,4][2]
p_grim_9res <- c(p_grim_9res, p_grim_9res_1)
coef_grim_9res_1 <- coef(summary(model_nn_9_grim))[,1][2]
coef_grim_9res <- c(coef_grim_9res, coef_grim_9res_1)

p_grim_res_1 <- coef(summary(model_nn_dif_grim))[,4][2]
p_grim_res <- c(p_grim_res, p_grim_res_1)
coef_grim_res_1 <- coef(summary(model_nn_dif_grim))[,1][2]
coef_grim_res <- c(coef_grim_res, coef_grim_res_1)

p_grim_pace_1 <- coef(summary(model_nn_pace_grim))[,4][2]
p_grim_pace <- c(p_grim_pace, p_grim_pace_1)
coef_grim_pace_1 <- coef(summary(model_nn_pace_grim))[,1][2]
coef_grim_pace <- c(coef_grim_pace, coef_grim_pace_1)

mean_grim_15res_exp <- c(mean_grim_15res_exp, mean(merged_both_nn_exp$age_15_grim_res))
mean_grim_15res_c <- c(mean_grim_15res_c, mean(merged_both_nn_c$age_15_grim_res))
mean_grim_res_exp <- c(mean_grim_res_exp, mean(merged_both_nn_exp$grim_dif))
mean_grim_res_c <- c(mean_grim_res_c, mean(merged_both_nn_c$grim_dif))
mean_grim_pace_exp <- c(mean_grim_pace_exp, mean(merged_both_nn_exp$grim_pace))
mean_grim_pace_c <- c(mean_grim_pace_c, mean(merged_both_nn_c$grim_pace))


p_dunedin_15res_1 <- coef(summary(model_nn_15_dunedin))[,4][2]
p_dunedin_15res <- c(p_dunedin_15res, p_dunedin_15res_1)
coef_dunedin_15res_1 <- coef(summary(model_nn_15_dunedin))[,1][2]
coef_dunedin_15res <- c(coef_dunedin_15res, coef_dunedin_15res_1)

p_dunedin_9res_1 <- coef(summary(model_nn_9_dunedin))[,4][2]
p_dunedin_9res <- c(p_dunedin_9res, p_dunedin_9res_1)
coef_dunedin_9res_1 <- coef(summary(model_nn_9_dunedin))[,1][2]
coef_dunedin_9res <- c(coef_dunedin_9res, coef_dunedin_9res_1)

p_dunedin_res_1 <- coef(summary(model_nn_dif_dunedin))[,4][2]
p_dunedin_res <- c(p_dunedin_res, p_dunedin_res_1)
coef_dunedin_res_1 <- coef(summary(model_nn_dif_dunedin))[,1][2]
coef_dunedin_res <- c(coef_dunedin_res, coef_dunedin_res_1)

p_dunedin_pace_1 <- coef(summary(model_nn_pace_dunedin))[,4][2]
p_dunedin_pace <- c(p_dunedin_pace, p_dunedin_pace_1)
coef_dunedin_pace_1 <- coef(summary(model_nn_pace_dunedin))[,1][2]
coef_dunedin_pace <- c(coef_dunedin_pace, coef_dunedin_pace_1)

mean_dunedin_15res_exp <- c(mean_dunedin_15res_exp, mean(merged_both_nn_exp$age_15_dunedin_res))
mean_dunedin_15res_c <- c(mean_dunedin_15res_c, mean(merged_both_nn_c$age_15_dunedin_res))
mean_dunedin_res_exp <- c(mean_dunedin_res_exp, mean(merged_both_nn_exp$dunedin_dif))
mean_dunedin_res_c <- c(mean_dunedin_res_c, mean(merged_both_nn_c$dunedin_dif))
mean_dunedin_pace_exp <- c(mean_dunedin_pace_exp, mean(merged_both_nn_exp$dunedin_pace))
mean_dunedin_pace_c <- c(mean_dunedin_pace_c, mean(merged_both_nn_c$dunedin_pace))




p_dnam_15res_1 <- coef(summary(model_nn_15_dnam))[,4][2]
p_dnam_15res <- c(p_dnam_15res, p_dnam_15res_1)
coef_dnam_15res_1 <- coef(summary(model_nn_15_dnam))[,1][2]
coef_dnam_15res <- c(coef_dnam_15res, coef_dnam_15res_1)

p_dnam_9res_1 <- coef(summary(model_nn_9_dnam))[,4][2]
p_dnam_9res <- c(p_dnam_9res, p_dnam_9res_1)
coef_dnam_9res_1 <- coef(summary(model_nn_9_dnam))[,1][2]
coef_dnam_9res <- c(coef_dnam_9res, coef_dnam_9res_1)

p_dnam_res_1 <- coef(summary(model_nn_dif_dnam))[,4][2]
p_dnam_res <- c(p_dnam_res, p_dnam_res_1)
coef_dnam_res_1 <- coef(summary(model_nn_dif_dnam))[,1][2]
coef_dnam_res <- c(coef_dnam_res, coef_dnam_res_1)

p_dnam_pace_1 <- coef(summary(model_nn_pace_dnam))[,4][2]
p_dnam_pace <- c(p_dnam_pace, p_dnam_pace_1)
coef_dnam_pace_1 <- coef(summary(model_nn_pace_dnam))[,1][2]
coef_dnam_pace <- c(coef_dnam_pace, coef_dnam_pace_1)

mean_dnam_15res_exp <- c(mean_dnam_15res_exp, mean(merged_both_nn_exp$age_15_dnam_res))
mean_dnam_15res_c <- c(mean_dnam_15res_c, mean(merged_both_nn_c$age_15_dnam_res))
mean_dnam_res_exp <- c(mean_dnam_res_exp, mean(merged_both_nn_exp$dnam_dif))
mean_dnam_res_c <- c(mean_dnam_res_c, mean(merged_both_nn_c$dnam_dif))
mean_dnam_pace_exp <- c(mean_dnam_pace_exp, mean(merged_both_nn_exp$dnam_pace))
mean_dnam_pace_c <- c(mean_dnam_pace_c, mean(merged_both_nn_c$dnam_pace))


p_ped_15res_1 <- coef(summary(model_nn_15_ped))[,4][2]
p_ped_15res <- c(p_ped_15res, p_ped_15res_1)
coef_ped_15res_1 <- coef(summary(model_nn_15_ped))[,1][2]
coef_ped_15res <- c(coef_ped_15res, coef_ped_15res_1)

p_ped_9res_1 <- coef(summary(model_nn_9_ped))[,4][2]
p_ped_9res <- c(p_ped_9res, p_ped_9res_1)
coef_ped_9res_1 <- coef(summary(model_nn_9_ped))[,1][2]
coef_ped_9res <- c(coef_ped_9res, coef_ped_9res_1)

p_ped_res_1 <- coef(summary(model_nn_dif_ped))[,4][2]
p_ped_res <- c(p_ped_res, p_ped_res_1)
coef_ped_res_1 <- coef(summary(model_nn_dif_ped))[,1][2]
coef_ped_res <- c(coef_ped_res, coef_ped_res_1)

p_ped_pace_1 <- coef(summary(model_nn_pace_ped))[,4][2]
p_ped_pace <- c(p_ped_pace, p_ped_pace_1)
coef_ped_pace_1 <- coef(summary(model_nn_pace_ped))[,1][2]
coef_ped_pace <- c(coef_ped_pace, coef_ped_pace_1)

p_hannum_15res_1 <- coef(summary(model_nn_15_hannum))[,4][2]
p_hannum_15res <- c(p_hannum_15res, p_hannum_15res_1)
coef_hannum_15res_1 <- coef(summary(model_nn_15_hannum))[,1][2]
coef_hannum_15res <- c(coef_hannum_15res, coef_hannum_15res_1)

p_hannum_9res_1 <- coef(summary(model_nn_9_hannum))[,4][2]
p_hannum_9res <- c(p_hannum_9res, p_hannum_9res_1)
coef_hannum_9res_1 <- coef(summary(model_nn_9_hannum))[,1][2]
coef_hannum_9res <- c(coef_hannum_9res, coef_hannum_9res_1)

p_hannum_res_1 <- coef(summary(model_nn_dif_hannum))[,4][2]
p_hannum_res <- c(p_hannum_res, p_hannum_res_1)
coef_hannum_res_1 <- coef(summary(model_nn_dif_hannum))[,1][2]
coef_hannum_res <- c(coef_hannum_res, coef_hannum_res_1)

p_hannum_pace_1 <- coef(summary(model_nn_pace_hannum))[,4][2]
p_hannum_pace <- c(p_hannum_pace, p_hannum_pace_1)
coef_hannum_pace_1 <- coef(summary(model_nn_pace_hannum))[,1][2]
coef_hannum_pace <- c(coef_hannum_pace, coef_hannum_pace_1)

mean_ped_15res_exp <- c(mean_ped_15res_exp, mean(merged_both_nn_exp$age_15_ped_res))
mean_ped_15res_c <- c(mean_ped_15res_c, mean(merged_both_nn_c$age_15_ped_res))
mean_ped_res_exp <- c(mean_ped_res_exp, mean(merged_both_nn_exp$ped_dif))
mean_ped_res_c <- c(mean_ped_res_c, mean(merged_both_nn_c$ped_dif))
mean_ped_pace_exp <- c(mean_ped_pace_exp, mean(merged_both_nn_exp$ped_pace))
mean_ped_pace_c <- c(mean_ped_pace_c, mean(merged_both_nn_c$ped_pace))


mean_hannum_15res_exp <- c(mean_hannum_15res_exp, mean(merged_both_nn_exp$age_15_hannum_res))
mean_hannum_15res_c <- c(mean_hannum_15res_c, mean(merged_both_nn_c$age_15_hannum_res))
mean_hannum_res_exp <- c(mean_hannum_res_exp, mean(merged_both_nn_exp$hannum_dif))
mean_hannum_res_c <- c(mean_hannum_res_c, mean(merged_both_nn_c$hannum_dif))
mean_hannum_pace_exp <- c(mean_hannum_pace_exp, mean(merged_both_nn_exp$hannum_pace))
mean_hannum_pace_c <- c(mean_hannum_pace_c, mean(merged_both_nn_c$hannum_pace))

p_bioage4_15res_1 <- coef(summary(model_nn_15_bioage4))[,4][2]
p_bioage4_15res <- c(p_bioage4_15res, p_bioage4_15res_1)
coef_bioage4_15res_1 <- coef(summary(model_nn_15_bioage4))[,1][2]
coef_bioage4_15res <- c(coef_bioage4_15res, coef_bioage4_15res_1)

p_bioage4_9res_1 <- coef(summary(model_nn_9_bioage4))[,4][2]
p_bioage4_9res <- c(p_bioage4_9res, p_bioage4_9res_1)
coef_bioage4_9res_1 <- coef(summary(model_nn_9_bioage4))[,1][2]
coef_bioage4_9res <- c(coef_bioage4_9res, coef_bioage4_9res_1)

p_bioage4_res_1 <- coef(summary(model_nn_dif_bioage4))[,4][2]
p_bioage4_res <- c(p_bioage4_res, p_bioage4_res_1)
coef_bioage4_res_1 <- coef(summary(model_nn_dif_bioage4))[,1][2]
coef_bioage4_res <- c(coef_bioage4_res, coef_bioage4_res_1)

p_bioage4_pace_1 <- coef(summary(model_nn_pace_bioage4))[,4][2]
p_bioage4_pace <- c(p_bioage4_pace, p_bioage4_pace_1)
coef_bioage4_pace_1 <- coef(summary(model_nn_pace_bioage4))[,1][2]
coef_bioage4_pace <- c(coef_bioage4_pace, coef_bioage4_pace_1)

mean_bioage4_15res_exp <- c(mean_bioage4_15res_exp, mean(merged_both_nn_exp$age_15_bioage4_res))
mean_bioage4_15res_c <- c(mean_bioage4_15res_c, mean(merged_both_nn_c$age_15_bioage4_res))
mean_bioage4_res_exp <- c(mean_bioage4_res_exp, mean(merged_both_nn_exp$bioage4_dif))
mean_bioage4_res_c <- c(mean_bioage4_res_c, mean(merged_both_nn_c$bioage4_dif))
mean_bioage4_pace_exp <- c(mean_bioage4_pace_exp, mean(merged_both_nn_exp$bioage4_pace))
mean_bioage4_pace_c <- c(mean_bioage4_pace_c, mean(merged_both_nn_c$bioage4_pace))


}

df_results <- data.frame(p_pack_15res, coef_pack_15res, p_pack_9res, coef_pack_9res, p_pack_res, coef_pack_res, p_pack_pace, coef_pack_pace, p_pheno_15res, coef_pheno_15res, p_pheno_9res, coef_pheno_9res, p_pheno_res, coef_pheno_res, p_pheno_pace, coef_pheno_pace, p_grim_15res, coef_grim_15res, p_grim_9res, coef_grim_9res, p_grim_res, coef_grim_res, p_grim_pace, coef_grim_pace, p_dunedin_15res, coef_dunedin_15res, p_dunedin_9res, coef_dunedin_9res, p_dunedin_res, coef_dunedin_res, p_dunedin_pace, coef_dunedin_pace, p_dnam_15res, coef_dnam_15res, p_dnam_9res, coef_dnam_9res, p_dnam_res, coef_dnam_res, p_dnam_pace, coef_dnam_pace, p_ped_15res, coef_ped_15res, p_ped_9res, coef_ped_9res, p_ped_res, coef_ped_res, p_ped_pace, coef_ped_pace, p_hannum_15res, coef_hannum_15res, p_hannum_9res, coef_hannum_9res, p_hannum_res, coef_hannum_res, p_hannum_pace, coef_hannum_pace, p_bioage4_15res, coef_bioage4_15res, p_bioage4_9res, coef_bioage4_res, p_bioage4_res, coef_bioage4_res, p_bioage4_pace, coef_bioage4_pace)

df_means <- data.frame(mean_pack_15res_exp, mean_pack_15res_c, mean_pack_res_exp, mean_pack_res_c, mean_pack_pace_exp, mean_pack_pace_c,mean_pheno_15res_exp, mean_pheno_15res_c, mean_pheno_res_exp, mean_pheno_res_c, mean_pheno_pace_exp, mean_pheno_pace_c, mean_grim_15res_exp, mean_grim_15res_c, mean_grim_res_exp, mean_grim_res_c, mean_grim_pace_exp, mean_grim_pace_c, mean_dunedin_15res_exp, mean_dunedin_15res_c, mean_dunedin_res_exp, mean_dunedin_res_c, mean_dunedin_pace_exp, mean_dunedin_pace_c, mean_dnam_15res_exp, mean_dnam_15res_c, mean_dnam_res_exp, mean_dnam_res_c, mean_dnam_pace_exp, mean_dnam_pace_c,mean_ped_15res_exp, mean_ped_15res_c, mean_ped_res_exp, mean_ped_res_c, mean_ped_pace_exp, mean_ped_pace_c, mean_hannum_15res_exp, mean_hannum_15res_c, mean_hannum_res_exp, mean_hannum_res_c, mean_hannum_pace_exp, mean_hannum_pace_c, mean_bioage4_15res_exp, mean_bioage4_15res_c, mean_bioage4_res_exp, mean_bioage4_res_c, mean_bioage4_pace_exp, mean_bioage4_pace_c)



write.csv(df_results, "/home/HPA/Hannah/clocks_slowmetab_alc_1.csv")
write.csv(df_means, "/home/HPA/Hannah/clocks_slowmetab_alc_1_means.csv")

