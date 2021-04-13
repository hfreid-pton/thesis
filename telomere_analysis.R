# telomere analysis
#code name on server: /home/HPA/Hannah/telomere_analysis1_correctTLfile_looped_1-27.R

library(dplyr)
library(tidyr)
library(stringr)

# making vectors to save pvals, coef, and means from 1000 iterations
p_val_newmethod_inc <- c()
coef_newmethod_inc <- c()
p_val_newmethod_dec <- c()
coef_newmethod_dec <- c()
p_val_newmethod_rate_inc <- c()
coef_newmethod_rate_inc <- c()
p_val_newmethod_rate_dec <- c()
coef_newmethod_rate_dec <- c()
p_val_reg15 <- c()
coef_reg15 <- c()
p_val_reg15_z <- c()
coef_reg15_z <- c()
p_val_reg9 <- c()
coef_reg9 <- c()
p_val_reg9_z <- c()
coef_reg9_z <- c()
p_val_reg_15_9 <- c()
coef_reg_15_9 <- c()
p_val_chi_dec <- c()
p_val_chi_inc <- c()
p_zscored_diff <- c()
coef_zscored_diff <- c()
p_diff_z <- c()
coef_diff_z <- c()
p_diff <- c()
coef_diff <- c()

mean_newmethodinc_exp <- c()
mean_newmethodinc_c <- c()
mean_newmethoddec_exp <- c()
mean_newmethoddec_c <- c()
mean_reg15_exp <- c()
mean_reg15_c <- c()
mean_reg9_exp <- c()
mean_reg9_c <- c()
mean_diff_exp <- c()
mean_diff_c <- c()
mean_diffz_exp <- c()
mean_diffz_c <- c()
prop_inc_c <- c()
prop_dec_c <- c()


# doing this 1000 times for nearest neighbor binary

# reading in survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

# reading in TL data
tl <- readRDS("/home/HPA/HPA_group_datasets/TLdata/TL_Sept2018.RDS")

tl <- tl %>% filter(Notes == "")

tl <- drop_na(tl)

# merging survey data with TL data
tl <- merge(tl, surv, by="ffFam")

# 1000 iterations (comparing exposed to 1000 iterations of matched controls)
for (i in 1:1000){
print(i)

#control matching (this is from MatchIt_11-1.R)

library(MatchIt)

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

#only keeping individuals in for matching that have tl data available
num_wells <- num_wells[num_wells$ffFam %in% tl$ffFam,]

exposed <- num_wells %>% filter(counts_10 > 0)


#finding potential controls to be used for matching
within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam

num_wells_10 <- num_wells
num_wells_10$exposure <- ifelse(num_wells_10$counts_10 > 0, 1, 0)
num_wells_10 <- num_wells_10[(num_wells_10$ffFam %in% controls_ffFam) | (num_wells_10$ffFam %in% exposed$ffFam),]

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

#keeping only individuals who didn't move or who moved and didn't change schools
merged4 <- merged4 %>% filter(p6j1 == 2 | p6j5 == 0)


merged4 <- merged4 %>% select(ffFam, cm1ethrace,exposure,cm1bsex, cm5povca, p5h15c, m1g4,k6d46, ck6pcgrel, cp6povca)

#changing cm1ethrace to dummy variables (white=reference)
merged4$his <- ifelse(merged4$cm1ethrace == 3, 1, 0)
merged4$afr <- ifelse(merged4$cm1ethrace == 2, 1, 0)

# matching
m.out <- matchit(exposure ~ his + afr + cm1bsex + cm5povca + p5h15c + m1g4 + k6d46 + ck6pcgrel + cp6povca, data=merged4, method="nearest", ratio=1)

match_df <- match.data(m.out)
intersect_ffFam <- intersect(tl$ffFam, match_df$ffFam)
match_df <- match_df[match_df$ffFam %in% intersect_ffFam,]
match_df_exp <- match_df %>% filter(exposure == 1)
match_df_control <- match_df %>% filter(exposure == 0)


# calculated rate of change in TL, difference in TL, etc.
tl$mean <- (tl$tTL + tl$meanadj_cTL)/2
tl$dif <- tl$tTL - tl$meanadj_cTL
tl <- tl %>% mutate(new_tl = (dif/mean))
tl <- tl %>% mutate(change_age = (ck6yagem - ch5agem))
tl <- tl %>% mutate(new_tl_rate = (new_tl/change_age))
tl <- tl %>% mutate(diff_z = (z_tTL - z_meanadjcTL))
tl <- tl %>% mutate(diff = tTL - (meanadj_cTL))
tl <- tl %>% mutate(zscored_diff = ((diff - mean(diff))/sd(diff)))


merged_tl <- merge(match_df, tl, by="ffFam")
#have to sep. increase and decrease in TL before taking the log because can't take log of neg number
newmethod_inc <- merged_tl %>% filter(new_tl > 0)
newmethod_dec <- merged_tl %>% filter(new_tl < 0)

reg_newmethod_inc <- glm(log(new_tl) ~ exposure + ck6yagem + ch5agem, data=newmethod_inc)
reg_newmethod_rate_inc <- glm(log(new_tl_rate) ~ exposure + ck6yagem + ch5agem, data=newmethod_inc)

reg_newmethod_dec <- glm(log(abs(new_tl)) ~ exposure + ck6yagem + ch5agem, data=newmethod_dec)
reg_newmethod_rate_dec <- glm(log(abs(new_tl_rate)) ~ exposure + ck6yagem + ch5agem, data=newmethod_dec)

reg_15 <- glm(log(tTL) ~ exposure + ck6yagem + ch5agem, data=merged_tl)
reg_9 <- glm(log(meanadj_cTL) ~ exposure + ck6yagem + ch5agem, data=merged_tl)

reg_15_z <- glm(z_tTL ~ exposure + ck6yagem + ch5agem, data=merged_tl)
reg_9_z <- glm(z_meanadjcTL ~ exposure + ck6yagem + ch5agem, data=merged_tl)

reg_15_9 <- glm(log(tTL) ~ exposure + log(meanadj_cTL) + ck6yagem + ch5agem, data=merged_tl)
reg_15_9_z <- glm(z_tTL ~ exposure + z_meanadjcTL + ck6yagem + ch5agem, data=merged_tl)

reg_diff <- glm(diff ~ exposure + ck6yagem + ch5agem, data=merged_tl)

reg_diff_z <- glm(diff_z ~ exposure + ck6yagem + ch5agem, data=merged_tl)
reg_zscored_diff <- glm(zscored_diff ~ exposure + ck6yagem + ch5agem, data=merged_tl)

# chi square for number of increased/decreased
diff_inc <- merged_tl %>% filter(diff > 0)
diff_dec <- merged_tl %>% filter(diff < 0)
diff_nochange <- merged_tl %>% filter(diff == 0)
exp_inc <- length(diff_inc$exposure[diff_inc$exposure == 1])
c_inc <- length(diff_inc$exposure[diff_inc$exposure == 0])
exp_dec <- length(diff_dec$exposure[diff_inc$exposure == 1])
c_dec <- length(diff_dec$exposure[diff_dec$exposure == 0])
exp_total <- length(merged_tl$exposure[merged_tl$exposure == 1])
c_total <- length(merged_tl$exposure[merged_tl$exposure == 0])
inc_chi <- prop.test(x=c(exp_inc, c_inc),n=c(exp_total, c_total))
dec_chi <- prop.test(x=c(exp_dec, c_dec), n=c(exp_total, c_total))


#saving pvals, coefficients, means from each iteration
p_newmethod_inc <- coef(summary(reg_newmethod_inc))[,4][2]
p_val_newmethod_inc <- c(p_val_newmethod_inc, p_newmethod_inc)
coef_newmethod1_inc <- coef(summary(reg_newmethod_inc))[,1][2]
coef_newmethod_inc <- c(coef_newmethod_inc, coef_newmethod1_inc)
p_newmethod_dec <- coef(summary(reg_newmethod_dec))[,4][2]
p_val_newmethod_dec <- c(p_val_newmethod_dec, p_newmethod_dec)
coef_newmethod1_dec <- coef(summary(reg_newmethod_dec))[,1][2]
coef_newmethod_dec <- c(coef_newmethod_dec, coef_newmethod1_dec)


p_newmethodrate_inc <- coef(summary(reg_newmethod_rate_inc))[,4][2]
p_val_newmethod_rate_inc <- c(p_val_newmethod_rate_inc, p_newmethodrate_inc)
coef_newmethod_rate1_inc <- coef(summary(reg_newmethod_rate_inc))[,1][2]
coef_newmethod_rate_inc <- c(coef_newmethod_rate_inc, coef_newmethod_rate1_inc)
p_newmethodrate_dec <- coef(summary(reg_newmethod_rate_dec))[,4][2]
p_val_newmethod_rate_dec <- c(p_val_newmethod_rate_dec, p_newmethodrate_dec)
coef_newmethod_rate1_dec <- coef(summary(reg_newmethod_rate_dec))[,1][2]
coef_newmethod_rate_dec <- c(coef_newmethod_rate_dec, coef_newmethod_rate1_dec)


p_reg15 <- coef(summary(reg_15))[,4][2]
p_val_reg15 <- c(p_val_reg15, p_reg15)
coef_reg15_1 <- coef(summary(reg_15))[,1][2]
coef_reg15 <- c(coef_reg15, coef_reg15_1)
p_reg15_z <- coef(summary(reg_15_z))[,4][2]
p_val_reg15_z <- c(p_val_reg15_z, p_reg15_z)
coef_reg15z_1 <- coef(summary(reg_15_z))[,1][2]
coef_reg15_z <- c(coef_reg15_z, coef_reg15z_1)
p_reg9 <- coef(summary(reg_9))[,4][2]
p_val_reg9 <- c(p_val_reg9, p_reg9)
coef_reg9_1 <- coef(summary(reg_9))[,1][2]
coef_reg9 <- c(coef_reg9, coef_reg9_1)
p_reg9_z <- coef(summary(reg_9_z))[,4][2]
p_val_reg9_z <- c(p_val_reg9_z, p_reg9_z)
coef_reg9z_1 <- coef(summary(reg_9_z))[,1][2]
coef_reg9_z <- c(coef_reg9_z, coef_reg9z_1)
p_reg15_9 <- coef(summary(reg_15_9))[,4][2]
p_val_reg_15_9 <- c(p_reg15_9, p_reg15_9)
coef_reg15_9_1 <- coef(summary(reg_15_9))[,1][2]
coef_reg_15_9 <- c(coef_reg_15_9, coef_reg15_9_1)

p_val_chi_dec <- c(p_val_chi_dec, dec_chi$p.value)
p_val_chi_inc <- c(p_val_chi_inc, inc_chi$p.value)

p_regdiff <- coef(summary(reg_diff))[,4][2]
p_diff <- c(p_diff, p_regdiff)
coef_regdiff <- coef(summary(reg_diff))[,1][2]
coef_diff <- c(coef_diff, coef_regdiff)

p_regdiffz <- coef(summary(reg_diff_z))[,4][2]
p_diff_z <- c(p_diff_z, p_regdiffz)
coef_regdiffz <- coef(summary(reg_diff_z))[,4][2]
coef_diff_z <- c(coef_diff_z, coef_regdiffz)

p_zscoreddiff <- coef(summary(reg_zscored_diff))[,4][2]
p_zscored_diff <- c(p_zscored_diff, p_zscoreddiff)
coef_zscoreddiff <- coef(summary(reg_zscored_diff))[,1][2]
coef_zscored_diff <- c(coef_zscored_diff, coef_zscoreddiff)

#means
newmethod_inc_exp <- newmethod_inc %>% filter(exposure == 1)
newmethod_inc_c <- newmethod_inc %>% filter(exposure == 0)
newmethod_dec_exp <- newmethod_dec %>% filter(exposure == 1)
newmethod_dec_c <- newmethod_dec %>% filter(exposure == 0)

mean_newmethodinc_exp <- c(mean_newmethodinc_exp, mean(newmethod_inc_exp$new_tl))
mean_newmethodinc_c <- c(mean_newmethodinc_c, mean(newmethod_inc_c$new_tl))
mean_newmethoddec_exp <- c(mean_newmethoddec_exp, mean(newmethod_dec_exp$new_tl))
mean_newmethoddec_c <- c(mean_newmethoddec_c, mean(newmethod_dec_c$new_tl))

merged_tl_exp <- merged_tl %>% filter(exposure == 1)
merged_tl_c <- merged_tl %>% filter(exposure == 0)

mean_reg15_exp <- c(mean_reg15_exp, mean(merged_tl_exp$tTL))
mean_reg15_c <- c(mean_reg15_c, mean(merged_tl_c$tTL))
mean_reg9_exp <- c(mean_reg9_exp, mean(merged_tl_exp$meanadj_cTL))
mean_reg9_c <- c(mean_reg9_c, mean(merged_tl_c$meanadj_cTL))
mean_diff_exp <- c(mean_diff_exp, mean(merged_tl_exp$diff))
mean_diff_c <- c(mean_diff_c, mean(merged_tl_c$diff))
mean_diffz_exp <- c(mean_diffz_exp, mean(merged_tl_exp$diff_z))
mean_diffz_c <- c(mean_diffz_c, mean(merged_tl_c$diff_z))



prop_inc_exp <- exp_inc/exp_total
prop_dec_exp <- exp_dec/exp_total
prop_inc_c <- c(prop_inc_c, (c_inc/c_total))
prop_dec_c <- c(prop_dec_c, (c_dec/c_total))

}



df_loop <- data.frame(p_val_newmethod_inc, coef_newmethod_inc, p_val_newmethod_dec, coef_newmethod_dec, p_val_newmethod_rate_inc, coef_newmethod_rate_inc, p_val_newmethod_rate_dec, coef_newmethod_rate_dec, p_val_reg15, coef_reg15, p_val_reg9, coef_reg9, p_val_reg_15_9, coef_reg_15_9, p_val_chi_dec, p_val_chi_inc, p_zscored_diff, coef_zscored_diff, p_diff_z, coef_diff_z, p_diff, coef_diff)

df_means <- data.frame(mean_newmethodinc_exp, mean_newmethodinc_c, mean_newmethoddec_exp, mean_newmethoddec_c, mean_reg15_exp, mean_reg15_c, mean_reg9_exp, mean_reg9_c, mean_diff_exp, mean_diff_c, mean_diffz_exp, mean_diffz_c, prop_inc_exp, prop_inc_c, prop_dec_exp, prop_dec_c)

write.csv(df_loop, "/home/HPA/Hannah/telomere_analysis_loop_3-21.csv")
write.csv(df_means, "/home/HPA/Hannah/telomere_analysis_loop_3-21_means.csv")
