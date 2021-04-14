# looking for association between health condition and metabolism rate in exposed group
# on server filename: health_metabolism_exposed_updated.R

library(dplyr)
library(tidyr)
library(stringr)


# setting up variables outside of loop to record p values and means (count variables are counts for k6d4 and bmi 
# but are proportions for other variables)
# see later in code for what each variable means
# exp=slower metab, c=faster metab (reusing original code, so kept notation)
# don't need a loop in this code because just run once (just one exposed group)
# but reusing code from original health scripts
p_k6d4 <- c()
count_k6d4_exp <- c()
count_k6d4_c <- c()

p_k6d8_earlier <- c()
count_k6d8_earlier_exp <- c()
count_k6d8_earlier_c <- c()

p_k6d8_later <- c()
count_k6d8_later_exp <- c()
count_k6d8_later_c <- c()

p_p6b5 <- c()
count_p6b5_exp <- c()
count_p6b5_c <- c()

p_p6b14 <- c()
count_p6b14_exp <- c()
count_p6b14_c <- c()

p_p6b8 <- c()
count_p6b8_exp <- c()
count_p6b8_c <- c()

p_p6b9_101 <- c()
count_p6b9_101_exp <- c()
count_p6b9_101_c <- c()

p_p6b19 <- c()
count_p6b19_exp <- c()
count_p6b19_c <- c()

p_p6b16 <- c()
count_p6b16_exp <- c()
count_p6b16_c <- c()

p_p6b23 <- c()
count_p6b23_exp <- c()
count_p6b23_c <- c()

p_asthma <- c()
count_asthma_exp <- c()
count_asthma_c <- c()

p_bmi <- c()
count_bmi_exp <- c()
count_bmi_c <- c()

p_bmiz <- c()
count_bmiz_exp <- c()
count_bmi_c <- c()


# reading in metabolism score data and creating binary scores
score <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/9var_scores_3-14.csv")

score$ffFam <- str_sub(score$Indiv, end=-5)
score$binary <- ifelse(score$total < 5,0,1)
score$total_alc <- (score$score_2e16 + score$score_adh1c + score$score_adh1b + score$score_e15)
score$binary_alc <- ifelse(score$total_alc < 2,0,1)

library(tidyr)
library(stringr)

#filtering scores to just fracking-exposed individuals
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")
dist_10 <- dist %>% filter(distance_ff_site <= 10)

ffFam_10 <- str_sub(unique(dist_10$ffid), end=-5)

score_exp <- intersect(score$ffFam, ffFam_10)

# reading in survey data
surv2 <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

#only keep people in surv that have scores
surv2 <- surv2[surv2$ffFam %in% score_exp,]

merged_score <- merge(score, surv2, by="ffFam")

#taking out people who moved and also moved schools
merged_score <- merged_score %>% filter(p6j1 == 2 | p6j5 == 0)


surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

# renamed this to continue using original code (see note below)
num_wells_10 <- merged_score

#to reuse code from original health analysis, I reuse exposure variable
# here exposure ==1 is a slower metabolizer, exposure==0 is a faster metabolizer
num_wells_10$exposure <- num_wells_10$binary #this line is changed to num_wells_10$binary_alc for alc score


#what is labeled exp later is slower metabolizers and c/control is faster 
# (using code from original analysis--harder to manipulate code on the server)

#days absent from school for health last yr
k6d4 <- c()

# development later or earlier than avg
k6d8 <- c()

#diagnosed with depression/anxiety
p6b5 <- c()

#eczema/skin allergies in last year
p6b14 <- c()

# seizures
p6b8 <- c()

#headaches/migraines
p6b9_101 <- c()

#trouble breathing
p6b19 <- c()

#frequent headaches past year
p6b16 <- c()


# seeing dr for illness
p6b23 <- c()

# asthma
p6b2 <- c()

# bmi
ch6bmiz <- c()
ck6cbmi <- c()


# added year 9 variables for bmi, asthma, diabetes
ch5bmiz <- c()
p5h1b <- c()
p5h3j <- c()


for (i in 1:length(num_wells_10$ffFam)){
merged2 <- merged[merged$ffFam == num_wells_10$ffFam[i], ]
k6d4_1 <- as.integer(merged2$k6d4[1])
k6d8_1 <- as.integer(merged2$k6d8[1])
p6b5_1 <- as.integer(merged2$p6b5[1])
p6b14_1 <- as.integer(merged2$p6b14[1])
p6b8_1 <- as.integer(merged2$p6b8[1])
p6b9_101_1 <- as.integer(merged2$p6b9_101[1])
p6b19_1 <- as.integer(merged2$p6b19[1])
p6b16_1 <- as.integer(merged2$p6b16[1])
p6b23_1 <- as.integer(merged2$p6b23[1])
p6b2_1 <- as.integer(merged2$p6b2[1])
ch6bmiz_1 <- as.integer(merged2$ch6bmiz[1])
ck6cbmi_1 <- as.integer(merged2$ck6cbmi[1])
ch5bmiz_1 <- as.integer(merged2$ch5bmiz[1])
p5h1b_1 <- as.integer(merged2$p5h1b[1])
p5h3j_1 <- as.integer(merged2$p5h3j[1])
k6d4 <- append(k6d4, k6d4_1)
k6d8 <- append(k6d8, k6d8_1)
p6b5 <- append(p6b5, p6b5_1)
p6b14 <- append(p6b14, p6b14_1)
p6b8 <- append(p6b8, p6b8_1)
p6b9_101 <- append(p6b9_101, p6b9_101)
p6b19 <- append(p6b19, p6b19_1)
p6b16 <- append(p6b16, p6b16_1)
p6b23 <- append(p6b23, p6b23_1)
p6b2 <- append(p6b2, p6b2_1)
ch6bmiz <- append(ch6bmiz, ch6bmiz_1)
ck6cbmi <- append(ck6cbmi, ck6cbmi_1)
ch5bmiz <- append(ch5bmiz, ch5bmiz_1)
p5h1b <- append(p5h1b, p5h1b_1)
p5h3j <- append(p5h3j, p5h3j_1)
}




num_wells_10$k6d4 <- k6d4
num_wells_10_k6d4 <- num_wells_10 %>% filter(k6d4 >= 0)
num_wells_10_k6d4_expo <- num_wells_10_k6d4 %>% filter(exposure == 1)
num_wells_10_k6d4_c <- num_wells_10_k6d4 %>% filter(exposure == 0)
k6d4_t2 <- t.test(num_wells_10_k6d4_expo$k6d4, num_wells_10_k6d4_c$k6d4, var.equal=TRUE)

p_k6d4 <- c(p_k6d4, k6d4_t2$p.value)
count_k6d4_exp <- c(count_k6d4_exp, mean(num_wells_10_k6d4_expo$k6d4))
count_k6d4_c <- c(count_k6d4_c, mean(num_wells_10_k6d4_c$k6d4))


num_wells_10$ck6cbmi <- ck6cbmi
num_wells_10_ck6cbmi <- num_wells_10 %>% filter(ck6cbmi >= 0)
num_wells_10_bmi_expo <- num_wells_10_ck6cbmi %>% filter(exposure == 1)
num_wells_10_bmi_c <- num_wells_10_ck6cbmi %>% filter(exposure == 0)
bmi_t2 <- t.test(num_wells_10_bmi_expo$ck6cbmi, num_wells_10_bmi_c$ck6cbmi, var.equal=TRUE)

p_bmi <- c(p_bmi, bmi_t2$p.value)
count_bmi_exp <- c(count_bmi_exp, mean(num_wells_10_bmi_expo$ck6cbmi))
count_bmi_c <- c(count_bmi_c, mean(num_wells_10_bmi_c$ck6cbmi))


num_wells_10$ch5bmiz <- ch5bmiz
num_wells_10$ch6bmiz <- ch6bmiz
num_wells_10_bmiz <- num_wells_10[is.na(num_wells_10$ch5bmiz) == F,]
num_wells_10_bmiz <- num_wells_10_bmiz[is.na(num_wells_10_bmiz$ch6bmiz) ==F,]
num_wells_10_bmiz$diff <- num_wells_10_bmiz$ch6bmiz - num_wells_10_bmiz$ch5bmiz
num_wells_10_bmiz_expo <- num_wells_10_bmiz %>% filter(exposure == 1)
num_wells_10_bmiz_c <- num_wells_10_bmiz %>% filter(exposure == 0)
bmiz_t2 <- t.test(num_wells_10_bmiz_expo$diff, num_wells_10_bmiz_c$diff, var.equal=TRUE)

p_bmiz <- c(p_bmiz, bmiz_t2$p.value)
count_bmiz_exp <- c(count_bmiz_exp, mean(num_wells_10_bmiz_expo$ck6cbmi))
count_bmiz_c <- c(count_bmiz_c, mean(num_wells_10_bmiz_c$ck6cbmi))


num_wells_10$k6d8 <- k6d8
# later dev
num_wells_10_k6d8_later <- num_wells_10 %>% filter(k6d8 >= 2)
num_wells_10_k6d8_later$k6d8_corrected <- ifelse(num_wells_10_k6d8_later$k6d8 == 2, 0, 1)

# earlier dev
num_wells_10_k6d8_earlier <- num_wells_10 %>% filter(k6d8 <= 2)
num_wells_10_k6d8_earlier <- num_wells_10 %>% filter(k6d8 >= 1)
num_wells_10_k6d8_earlier$k6d8_corrected <- ifelse(num_wells_10_k6d8_earlier$k6d8 == 2, 0, 1)

num_wells_10_k6d8_later_expo <- num_wells_10_k6d8_later %>% filter(exposure == 1)
num_wells_10_k6d8_later_c <- num_wells_10_k6d8_later %>% filter(exposure == 0)
num_wells_10_k6d8_earlier_expo <- num_wells_10_k6d8_earlier %>% filter(exposure == 1)
num_wells_10_k6d8_earlier_c <- num_wells_10_k6d8_earlier %>% filter(exposure== 0)

#for later
num_wells_10_k6d8_later_expo <- num_wells_10_k6d8_later_expo$k6d8_corrected[is.na(num_wells_10_k6d8_later_expo$k6d8_corrected) == F]
num_wells_10_k6d8_later_c <- num_wells_10_k6d8_later_c$k6d8_corrected[is.na(num_wells_10_k6d8_later_c$k6d8_corrected) == F]
exp_k6d8_later_10 <- c(length(num_wells_10_k6d8_later_expo[num_wells_10_k6d8_later_expo == 1]), length(num_wells_10_k6d8_later_expo[num_wells_10_k6d8_later_expo == 0]))
c_k6d8_later_10 <- c(length(num_wells_10_k6d8_later_c[num_wells_10_k6d8_later_c == 1]), length(num_wells_10_k6d8_later_c[num_wells_10_k6d8_later_c == 0]))

counts_k6d8_later_corrected <- fisher.test(as.table(rbind(exp_k6d8_later_10, c_k6d8_later_10)))



num_wells_10_k6d8_earlier_expo <- num_wells_10_k6d8_earlier_expo$k6d8_corrected[is.na(num_wells_10_k6d8_earlier_expo$k6d8_corrected) == F]
num_wells_10_k6d8_earlier_c <- num_wells_10_k6d8_earlier_c$k6d8_corrected[is.na(num_wells_10_k6d8_earlier_c$k6d8_corrected) == F]
exp_k6d8_earlier_10 <- c(length(num_wells_10_k6d8_earlier_expo[num_wells_10_k6d8_earlier_expo == 1]), length(num_wells_10_k6d8_earlier_expo[num_wells_10_k6d8_earlier_expo == 0]))
c_k6d8_earlier_10 <- c(length(num_wells_10_k6d8_earlier_c[num_wells_10_k6d8_earlier_c == 1]), length(num_wells_10_k6d8_earlier_c[num_wells_10_k6d8_earlier_c == 0]))

counts_k6d8_earlier_corrected <- fisher.test(as.table(rbind(exp_k6d8_earlier_10, c_k6d8_earlier_10)))



p_k6d8_earlier <- c(p_k6d8_earlier, counts_k6d8_earlier_corrected$p.value)
p_k6d8_later <- c(p_k6d8_later, counts_k6d8_later_corrected$p.value)


k6d8_later_c <- (c_k6d8_later_10[1]/sum(c_k6d8_later_10))
k6d8_later_exp <- (exp_k6d8_later_10[1]/sum(exp_k6d8_later_10))
k6d8_earlier_c <- (c_k6d8_earlier_10[1]/sum(c_k6d8_earlier_10))
k6d8_earlier_exp <- (exp_k6d8_earlier_10[1]/sum(exp_k6d8_earlier_10))

count_k6d8_later_c <- c(count_k6d8_later_c, k6d8_later_c)
count_k6d8_later_exp <- c(count_k6d8_later_exp, k6d8_later_exp)
count_k6d8_earlier_c <- c(count_k6d8_earlier_c, k6d8_earlier_c)
count_k6d8_earlier_exp <- c(count_k6d8_earlier_exp, k6d8_earlier_exp)




# add p6b5 values
num_wells_10$p6b5 <- p6b5

# get rid of all error codes
num_wells_10_p6b5 <- num_wells_10 %>% filter(p6b5 >= 0)

# make into 0,1 binary
num_wells_10_p6b5$p6b5_corrected <- ifelse(num_wells_10_p6b5$p6b5 == 2, 0, 1)

# df for control(faster metab) vs. exposed (slower metab)
num_wells_10_p6b5_expo <- num_wells_10_p6b5 %>% filter(exposure == 1)
num_wells_10_p6b5_c <- num_wells_10_p6b5 %>% filter(exposure == 0)

# removing NAs
num_wells_10_p6b5_expo <- num_wells_10_p6b5_expo$p6b5_corrected[is.na(num_wells_10_p6b5_expo$p6b5_corrected) == F]
num_wells_10_p6b5_c <- num_wells_10_p6b5_c$p6b5_corrected[is.na(num_wells_10_p6b5_c$p6b5_corrected) == F]

#creating vectors of with condition vs. without for both slow and fast
exp_p6b5_10 <- c(length(num_wells_10_p6b5_expo[num_wells_10_p6b5_expo == 1]), length(num_wells_10_p6b5_expo[num_wells_10_p6b5_expo == 0]))
c_p6b5_10 <- c(length(num_wells_10_p6b5_c[num_wells_10_p6b5_c == 1]), length(num_wells_10_p6b5_c[num_wells_10_p6b5_c == 0]))


counts_p6b5_corrected <- fisher.test(as.table(rbind(exp_p6b5_10, c_p6b5_10)))


p_p6b5 <- c(p_p6b5, counts_p6b5_corrected$p.value)


p6b5_c <- (c_p6b5_10[1]/sum(c_p6b5_10))
p6b5_exp <- (exp_p6b5_10[1]/sum(exp_p6b5_10))

count_p6b5_c <- c(count_p6b5_c, p6b5_c)
count_p6b5_exp <- c(count_p6b5_exp, p6b5_exp)



num_wells_10$p6b14 <- p6b14
num_wells_10_p6b14 <- num_wells_10 %>% filter(p6b14 >= 0)
num_wells_10_p6b14$p6b14_corrected <- ifelse(num_wells_10_p6b14$p6b14 == 2, 0, 1)

num_wells_10_p6b14_expo <- num_wells_10_p6b14 %>% filter(exposure == 1)
num_wells_10_p6b14_c <- num_wells_10_p6b14 %>% filter(exposure == 0)
num_wells_10_p6b14_expo <- num_wells_10_p6b14_expo$p6b14_corrected[is.na(num_wells_10_p6b14_expo$p6b14_corrected) == F]
num_wells_10_p6b14_c <- num_wells_10_p6b14_c$p6b14_corrected[is.na(num_wells_10_p6b14_c$p6b14_corrected) == F]
exp_p6b14_10 <- c(length(num_wells_10_p6b14_expo[num_wells_10_p6b14_expo == 1]), length(num_wells_10_p6b14_expo[num_wells_10_p6b14_expo == 0]))
c_p6b14_10 <- c(length(num_wells_10_p6b14_c[num_wells_10_p6b14_c == 1]), length(num_wells_10_p6b14_c[num_wells_10_p6b14_c == 0]))

counts_p6b14_corrected <- fisher.test(as.table(rbind(exp_p6b14_10, c_p6b14_10)))


prop_c <- (c_p6b14_10[1]/sum(c_p6b14_10))
prop_c <- prop_c * sum(exp_p6b14_10)
prop_c_p6b14 <- c(prop_c, (sum(exp_p6b14_10)-prop_c))


p_p6b14 <- c(p_p6b14, counts_p6b14_corrected$p.value)

p6b14_c <- (c_p6b14_10[1]/sum(c_p6b14_10))
p6b14_exp <- (exp_p6b14_10[1]/sum(exp_p6b14_10))

count_p6b14_c <- c(count_p6b14_c, p6b14_c)
count_p6b14_exp <- c(count_p6b14_exp, p6b14_exp)




num_wells_10$p6b8 <- p6b8
num_wells_10_p6b8 <- num_wells_10 %>% filter(p6b8 >= 1)
num_wells_10_p6b8$p6b8_corrected <- ifelse(num_wells_10_p6b8$p6b8 == 2, 0, 1)

num_wells_10_p6b8_expo <- num_wells_10_p6b8 %>% filter(exposure == 1)
num_wells_10_p6b8_c <- num_wells_10_p6b8 %>% filter(exposure == 0)
num_wells_10_p6b8_expo <- num_wells_10_p6b8_expo$p6b8_corrected[is.na(num_wells_10_p6b8_expo$p6b8_corrected) == F]
num_wells_10_p6b8_c <- num_wells_10_p6b8_c$p6b8_corrected[is.na(num_wells_10_p6b8_c$p6b8_corrected) == F]
exp_p6b8_10 <- c(length(num_wells_10_p6b8_expo[num_wells_10_p6b8_expo == 1]), length(num_wells_10_p6b8_expo[num_wells_10_p6b8_expo == 0]))
c_p6b8_10 <- c(length(num_wells_10_p6b8_c[num_wells_10_p6b8_c == 1]), length(num_wells_10_p6b8_c[num_wells_10_p6b8_c == 0]))

counts_p6b8_corrected <- fisher.test(as.table(rbind(exp_p6b8_10, c_p6b8_10)))

p_p6b8 <- c(p_p6b8, counts_p6b8_corrected$p.value)

p6b8_c <- (c_p6b8_10[1]/sum(c_p6b8_10))
p6b8_exp <- (exp_p6b8_10[1]/sum(exp_p6b8_10))

count_p6b8_c <- c(count_p6b8_c, p6b8_c)
count_p6b8_exp <- c(count_p6b8_exp, p6b8_exp)


num_wells_10$p6b19 <- p6b19
num_wells_10_p6b19 <- num_wells_10 %>% filter(p6b19 >= 1)
num_wells_10_p6b19$p6b19_corrected <- ifelse(num_wells_10_p6b19$p6b19 == 2, 0, 1)

num_wells_10_p6b19_expo <- num_wells_10_p6b19 %>% filter(exposure == 1)
num_wells_10_p6b19_c <- num_wells_10_p6b19 %>% filter(exposure == 0)
num_wells_10_p6b19_expo <- num_wells_10_p6b19_expo$p6b19_corrected[is.na(num_wells_10_p6b19_expo$p6b19_corrected) == F]
num_wells_10_p6b19_c <- num_wells_10_p6b19_c$p6b19_corrected[is.na(num_wells_10_p6b19_c$p6b19_corrected) == F]
exp_p6b19_10 <- c(length(num_wells_10_p6b19_expo[num_wells_10_p6b19_expo == 1]), length(num_wells_10_p6b19_expo[num_wells_10_p6b19_expo == 0]))
c_p6b19_10 <- c(length(num_wells_10_p6b19_c[num_wells_10_p6b19_c == 1]), length(num_wells_10_p6b19_c[num_wells_10_p6b19_c == 0]))

counts_p6b19_corrected <- fisher.test(as.table(rbind(exp_p6b19_10, c_p6b19_10)))


p_p6b19 <- c(p_p6b19, counts_p6b19_corrected$p.value)

p6b19_c <- (c_p6b19_10[1]/sum(c_p6b19_10))
p6b19_exp <- (exp_p6b19_10[1]/sum(exp_p6b19_10))

count_p6b19_c <- c(count_p6b19_c, p6b19_c)
count_p6b19_exp <- c(count_p6b19_exp, p6b19_exp)



num_wells_10$p6b16 <- p6b16
num_wells_10_p6b16 <- num_wells_10 %>% filter(p6b16 >= 1)
num_wells_10_p6b16$p6b16_corrected <- ifelse(num_wells_10_p6b16$p6b16 == 2,0,1)
num_wells_10_p6b16_expo <- num_wells_10_p6b16 %>% filter(exposure == 1)
num_wells_10_p6b16_c <- num_wells_10_p6b16 %>% filter(exposure == 0)
num_wells_10_p6b16_expo <- num_wells_10_p6b16_expo$p6b16_corrected[is.na(num_wells_10_p6b16_expo$p6b16_corrected) == F]
num_wells_10_p6b16_c <- num_wells_10_p6b16_c$p6b16_corrected[is.na(num_wells_10_p6b16_c$p6b16_corrected) == F]
exp_p6b16_10 <- c(length(num_wells_10_p6b16_expo[num_wells_10_p6b16_expo == 1]), length(num_wells_10_p6b16_expo[num_wells_10_p6b16_expo == 0]))
c_p6b16_10 <- c(length(num_wells_10_p6b16_c[num_wells_10_p6b16_c == 1]), length(num_wells_10_p6b16_c[num_wells_10_p6b16_c == 0]))

counts_p6b16_corrected <- fisher.test(as.table(rbind(exp_p6b16_10, c_p6b16_10)))

p_p6b16 <- c(p_p6b16, counts_p6b16_corrected$p.value)


p6b16_c <- (c_p6b16_10[1]/sum(c_p6b16_10))
p6b16_exp <- (exp_p6b16_10[1]/sum(exp_p6b16_10))

count_p6b16_c <- c(count_p6b16_c, p6b16_c)
count_p6b16_exp <- c(count_p6b16_exp, p6b16_exp)


num_wells_10$p6b23 <- p6b23
num_wells_10_p6b23 <- num_wells_10 %>% filter(p6b23 >= 1)
num_wells_10_p6b23$p6b23_corrected <- ifelse(num_wells_10_p6b23$p6b23 == 2, 0, 1)
num_wells_10_p6b23_expo <- num_wells_10_p6b23 %>% filter(exposure == 1)
num_wells_10_p6b23_c <- num_wells_10_p6b23 %>% filter(exposure == 0)
num_wells_10_p6b23_expo <- num_wells_10_p6b23_expo$p6b23_corrected[is.na(num_wells_10_p6b23_expo$p6b23_corrected) == F]
num_wells_10_p6b23_c <- num_wells_10_p6b23_c$p6b23_corrected[is.na(num_wells_10_p6b23_c$p6b23_corrected) == F]
exp_p6b23_10 <- c(length(num_wells_10_p6b23_expo[num_wells_10_p6b23_expo == 1]), length(num_wells_10_p6b23_expo[num_wells_10_p6b23_expo == 0]))
c_p6b23_10 <- c(length(num_wells_10_p6b23_c[num_wells_10_p6b23_c == 1]), length(num_wells_10_p6b23_c[num_wells_10_p6b23_c == 0]))

counts_p6b23_corrected <- fisher.test(as.table(rbind(exp_p6b23_10, c_p6b23_10)))

p_p6b23 <- c(p_p6b23, counts_p6b23_corrected$p.value)

p6b23_c <- (c_p6b23_10[1]/sum(c_p6b23_10))
p6b23_exp <- (exp_p6b23_10[1]/sum(exp_p6b23_10))

count_p6b23_c <- c(count_p6b23_c, p6b23_c)
count_p6b23_exp <- c(count_p6b23_exp, p6b23_exp)



num_wells_10$p6b2 <- p6b2
num_wells_10$p5h1b <- p5h1b
num_wells_10_p6b2 <- num_wells_10 %>% filter((p6b2 >= 1) & (p5h1b >= 1))
num_wells_10_p6b2$p6b2_corrected <- ifelse(num_wells_10_p6b2$p6b2 == 2, 0, 1)
num_wells_10_p6b2$p5h1b_corrected <- ifelse(num_wells_10_p6b2$p5h1b == 2,0,1)
num_wells_10_p6b2$difference <- num_wells_10_p6b2$p6b2_corrected - num_wells_10_p6b2$p5h1b_corrected
num_wells_10_p6b2_expo <- num_wells_10_p6b2 %>% filter(exposure == 1)
num_wells_10_p6b2_c <- num_wells_10_p6b2 %>% filter(exposure == 0)
num_wells_10_p6b2_expo <- num_wells_10_p6b2_expo$difference
num_wells_10_p6b2_c <- num_wells_10_p6b2_c$difference
exp_p6b2_10 <- c(length(num_wells_10_p6b2_expo[num_wells_10_p6b2_expo == 1]), length(num_wells_10_p6b2_expo[num_wells_10_p6b5_expo == 0]))
c_p6b2_10 <- c(length(num_wells_10_p6b2_c[num_wells_10_p6b2_c == 1]), length(num_wells_10_p6b2_c[num_wells_10_p6b2_c == 0]))

counts_p6b2_corrected <- fisher.test(as.table(rbind(exp_p6b2_10, c_p6b2_10)))


p_asthma <- c(p_asthma, counts_p6b2_corrected$p.value)

p6b2_c <- (c_p6b2_10[1]/sum(c_p6b2_10))
p6b2_exp <- (exp_p6b2_10[1]/sum(exp_p6b2_10))

count_asthma_c <- c(count_asthma_c, p6b2_c)
count_asthma_exp <- c(count_asthma_exp, p6b2_exp)



num_wells_10$p6b6 <- p6b6
num_wells_10$p5h3j <- p5h3j
num_wells_10_p6b6 <- num_wells_10 %>% filter((p6b6 >= 1) & (p5h3j >= 1))
num_wells_10_p6b6$p6b6_corrected <- ifelse(num_wells_10_p6b6$p6b6 == 2, 0, 1)
num_wells_10_p6b6$p5h3j_corrected <- ifelse(num_wells_10_p6b6$p5h3j == 2, 0, 1)
num_wells_10_p6b6$difference <- num_wells_10_p6b6$p6b6_corrected - num_wells_10_p6b6$p5h3j_corrected
num_wells_10_p6b6_expo <- num_wells_10_p6b6 %>% filter(exposure == 1)
num_wells_10_p6b6_c <- num_wells_10_p6b6 %>% filter(exposure == 0)
num_wells_10_p6b6_expo <- num_wells_10_p6b6_expo$difference
num_wells_10_p6b6_c <- num_wells_10_p6b6_c$difference
exp_p6b6_10 <- c(length(num_wells_10_p6b6_expo[num_wells_10_p6b6_expo == 1]), length(num_wells_10_p6b6_expo[num_wells_10_p6b6_expo == 0]))
c_p6b6_10 <- c(length(num_wells_10_p6b6_c[num_wells_10_p6b6_c == 1]), length(num_wells_10_p6b6_c[num_wells_10_p6b6_c == 0]))

counts_p6b6_corrected <- fisher.test(as.table(rbind(exp_p6b6_10, c_p6b6_10)))

}

df_loop <- data.frame(p_k6d8_earlier, count_k6d8_earlier_exp, count_k6d8_earlier_c, p_k6d8_later, count_k6d8_later_exp, count_k6d8_later_c, p_p6b5, count_p6b5_exp, count_p6b5_c, p_p6b14, count_p6b14_exp, count_p6b14_c, p_p6b8, count_p6b8_exp, count_p6b8_c, p_p6b19, count_p6b19_exp, count_p6b19_c, p_p6b16, count_p6b16_exp, count_p6b16_c, p_p6b23, count_p6b23_exp, count_p6b23_c, p_asthma, count_asthma_exp, count_asthma_c)

write.csv(df_loop, "/home/HPA/Hannah/health_metabolism_exposed_all_1.csv")


