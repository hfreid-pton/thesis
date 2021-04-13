#health analysis--bmi z-score
#server code name: /home/HPA/Hannah/chisq_10_newmatch_nomove_loop_bmi2.R

library(dplyr)
library(tidyr)



p_bmi <- c()

p_bmi2 <- c()

p_bmidif <- c()

count_bmi_exp <- c()
count_bmi_c <- c()
count_bmidif_exp <- c()
count_bmidif_c <- c()
count_bmi2_exp <- c()
count_bmi2_c <- c()


# 1000 iterations of control matching

for (n in 1:1000){
print(n)

#matching
source("/home/HPA/Hannah/matchit_11-1.R")
match_df <- match.data(m.out)
match_df_exp <- match_df %>% filter(exposure == 1)
match_df_control <- match_df %>% filter(exposure == 0)

exposed_ffFam <- match_df_exp$ffFam
control_ffFam <- match_df_control$ffFam

#reading in survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

num_wells_10 <- match_df


ch5bmiz <- c()
ch6bmiz <- c()

#making df of bmi variables for exposed and matched controls
for (i in 1:length(num_wells_10$ffFam)){
merged2 <- merged[merged$ffFam == num_wells_10$ffFam[i], ]
ch6bmiz_1 <- as.integer(merged2$ch6bmiz[1])
ch5bmiz_1 <- as.integer(merged2$ch5bmiz[1])
ch5bmiz <- append(ch5bmiz, ch5bmiz_1)
ch6bmiz <- append(ch6bmiz, ch6bmiz_1)
}


#comparing year 15 zscored bmi
num_wells_10$ch6bmiz <- ch6bmiz
num_wells_10_bmi_expo <- num_wells_10 %>% filter(exposure == 1)
num_wells_10_bmi_c <- num_wells_10 %>% filter(exposure == 0)
num_wells_10_bmi_expo <- num_wells_10_bmi_expo$ch6bmiz[is.na(num_wells_10_bmi_expo$ch6bmiz) == F]
num_wells_10_bmi_c <- num_wells_10_bmi_c$ch6bmiz[is.na(num_wells_10_bmi_c$ch6bmiz) == F]
bmi_t <- t.test(num_wells_10_bmi_expo, num_wells_10_bmi_c)
bmi_t2 <- t.test(num_wells_10_bmi_expo, num_wells_10_bmi_c, var.equal=TRUE)



p_bmi <- c(p_bmi, bmi_t2$p.value)


count_bmi_exp <- c(count_bmi_exp, mean(num_wells_10_bmi_expo))
count_bmi_c <- c(count_bmi_c, mean(num_wells_10_bmi_c))

#comparing year 9 zscored bmi
num_wells_10$ch5bmiz <- ch5bmiz
num_wells_10_bmi2_expo <- num_wells_10 %>% filter(exposure == 1)
num_wells_10_bmi2_c <- num_wells_10 %>% filter(exposure == 0)
num_wells_10_bmi2_expo <- num_wells_10_bmi2_expo$ch5bmiz[is.na(num_wells_10_bmi2_expo$ch5bmiz) == F]
num_wells_10_bmi2_c <- num_wells_10_bmi2_c$ch5bmiz[is.na(num_wells_10_bmi2_c$ch5bmiz) == F]
bmi2_t <- t.test(num_wells_10_bmi2_expo, num_wells_10_bmi2_c)
bmi2_t2 <- t.test(num_wells_10_bmi2_expo, num_wells_10_bmi2_c, var.equal=TRUE)

p_bmi2 <- c(p_bmi2, bmi2_t2$p.value)

count_bmi2_exp <- c(count_bmi2_exp, mean(num_wells_10_bmi2_expo))
count_bmi2_c <- c(count_bmi2_c, mean(num_wells_10_bmi2_c))

#comparing change in zscored bmi
num_wells_10_dif <- num_wells_10[is.na(num_wells_10$ch5bmiz) == F,]
num_wells_10_dif <- num_wells_10_dif[is.na(num_wells_10_dif$ch6bmiz) == F,]
num_wells_10_dif$dif <- num_wells_10_dif$ch6bmiz - num_wells_10_dif$ch5bmiz
num_wells_10_dif_expo <- num_wells_10_dif %>% filter(exposure == 1)
num_wells_10_dif_c <- num_wells_10_dif %>% filter(exposure == 0)

bmidif_t2 <- t.test(num_wells_10_dif_expo$dif, num_wells_10_dif_c$dif, var.equal=TRUE)

p_bmidif <- c(p_bmidif, bmidif_t2$p.value)

count_bmidif_exp <- c(count_bmidif_exp, mean(num_wells_10_dif_expo$dif))
count_bmidif_c <- c(count_bmidif_c, mean(num_wells_10_dif_c$dif))



}

df_bmi <- data.frame(p_bmi,count_bmi_exp, count_bmi_c, p_bmi2, count_bmi2_exp, count_bmi2_c,p_bmidif, count_bmidif_exp, count_bmidif_c)
write.csv(df_bmi, "/home/HPA/Hannah/bmi_z_loop_results_3-9.csv")
