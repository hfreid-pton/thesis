##Many portions/much of this code is based on and adapted from code written by Ben Jacobson '19
# for his senior thesis
# See citation for his thesis in my written thesis References section

# finding sites (on 450k array) that are differentially methylated from age 9 to age 15
# in exposed group and separately in 100 iterations of matched control groups
# code file name on server: /home/HPA/Hannah/methylation_1_test_450k_2.R

library(dplyr)
library(minfi)
library(tidyr)
library(missMethyl)
library(gtools)
library(limma)
library(MatchIt)
library(stringr)

#450k Data
samp <- readRDS("/home/HPA/HPA_group_datasets/MethylationData/April2019/k450_April2019_preprocessNOOB_BMIQ/450K_preprocessNOOB_BMIQsamplesheet.RDS")
prob <- readRDS("/home/HPA/Hannah/BenCode/datafiles/probIDs450.RDS")
beta <- readRDS("/home/HPA/HPA_group_datasets/MethylationData/April2019/k450_April2019_preprocessNOOB_BMIQ/450K_all_preprocessnoobbmiqIDfixed.RDS")

#added problem samples
newprob <- c("F641205","F652256","F827700", "F842982", "F863570", "F925575", "F970208")
prob <- c(prob, newprob)
prob <- unique(prob)


# IDs started with a number, which caused an error
# added X to start of MethID in samp to fix error
samp$MethID <- paste("X", samp$MethID, sep="")


# remove SNP cpgs
Zhou_cpg <- readRDS("/home/HPA/HPA_group_datasets/MethylationData/ProbeLists/Zhou/hm450.hg19.manifest.rds")
Zhou_cpg_2 <- as(Zhou_cpg, "data.frame")
Zhou_cpg_2$cpg <- rownames(Zhou_cpg_2)
cpg_keep <- Zhou_cpg_2[!Zhou_cpg_2$MASK_general,]
cpg_keep_list <- cpg_keep$cpg

beta <- beta[rownames(beta) %in% cpg_keep_list, ]



#winsorizing
beta[beta > 0.99] <- 0.99
beta[beta < 0.01] <- 0.01

#remove problem samples
samp <- samp[!(samp$ffFam %in% prob),]

samp_9 <- samp %>% filter(Group == "FC09")
samp_15 <- samp %>% filter(Group == "FC15")

# selecting only individuals with samples at 9 and 15
samp_both_ffFam <- intersect(samp_9$ffFam, samp_15$ffFam)

samp <- samp[samp$ffFam %in% samp_both_ffFam,]

# read INC files
INC <- readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p1.RDS")
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p2.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p3.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p4.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p5.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p6.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p7.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p8.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p9.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p10.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p11.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p12.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p13.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p14.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p15.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p16.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p17.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p18.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p19.RDS"))
INC <- cbind(INC, readRDS("/home/HPA/Jacobson/Lisa/k450_April2019_preprocessNOOB_BMIQ/INCs_450K/INCs_p23.RDS"))


#fixing error discussed above by adding X
colnames(INC) <- paste("X", colnames(INC), sep="")

# set up vectors outside of loop to store data from each iteration of control EWAS
cpg <- c()
cpg_9 <- c()
cpg_c <- c()
F.p <- c()
F.p_9 <- c()
F.p_c <- c()
F.p.BH <- c()
F.p.BH_9 <- c()
F.p.BH_c <- c()
p_X1.1 <- c()
p_X1.1_9 <- c()
p_X1.1_c <- c()
p.BH_X1.1 <- c()
p.BH_X1.1_9 <- c()
p.BH_X1.1_c <- c()
b_X1.1 <- c()
b_X1.1_9 <- c()
b_X1.1_c <- c()
sigma2 <- c()
sigma2_9 <- c()
sigma2_c <- c()
var.b_X1.1 <- c()
var.b_X1.1_9 <- c()
var.b_X1.1_c <- c()
fit.ctl <- c()
fit.ctl_9 <- c()
fit.ctl_c <- c()
mean <- c()
mean_9 <- c()
mean_c <- c()


#reading in fracking distance file (for matching)
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")
#reading in survey data (for matching)
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")


for (n in 1:100){
print(n)

# matching controls
#defining exposed individuals
dist_10 <- dist %>% filter(distance_ff_site <= 10)
all_ffid <- unique(dist$ffid)


ffid_10 <- dist_10$ffid
counts_10 <- c()

for (m in 1:length(all_ffid)){
num_10 <- length(ffid_10[ffid_10 == all_ffid[m]])
counts_10 <- append(counts_10, num_10)
}

num_wells <- data.frame(all_ffid, counts_10)

#changing ffid to ffFam id
ffam_test <- str_sub(num_wells$all_ffid, end=-5)
num_wells$ffFam <- ffam_test

#removing individuals not in 450k data
num_wells <- num_wells[num_wells$ffFam %in% samp$ffFam,]

exposed <- num_wells %>% filter(counts_10 > 0)

# defining potential controls for matching
within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam

num_wells_10 <- num_wells
num_wells_10$exposure <- ifelse(num_wells_10$counts_10 > 0, 1, 0)
num_wells_10 <- num_wells_10[(num_wells_10$ffFam %in% controls_ffFam) | (num_wells_10$ffFam %in% exposed$ffFam),]

merged2 <- surv[surv$ffFam %in% num_wells_10$ffFam,]

merged3 <- merge(merged2, num_wells_10, by="ffFam")

# getting rid of samples with NAs or error codes for variables used in matching
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

#keeping people who didn't move or who moved but didn't switch schools
merged4 <- merged4 %>% filter(p6j1 == 2 | p6j5 == 0)


merged4 <- merged4 %>% select(ffFam, cm1ethrace,exposure,cm1bsex, cm5povca, p5h15c, m1g4,k6d46, ck6pcgrel, cp6povca)

#making cm1ethrace into dummary var
merged4$his <- ifelse(merged4$cm1ethrace == 3, 1, 0)
merged4$afr <- ifelse(merged4$cm1ethrace == 2, 1, 0)

#matching
m.out <- matchit(exposure ~ his + afr + cm1bsex + cm5povca + p5h15c + m1g4 + k6d46 + ck6pcgrel + cp6povca, data=merged4, method="nearest", ratio=1)


match_df <- match.data(m.out)
match_df_exp <- match_df %>% filter(exposure == 1)
match_df_control <- match_df %>% filter(exposure == 0)

exposed_ffFam <- match_df_exp$ffFam
control_ffFam <- match_df_control$ffFam


#selecting only people in matchdf
samp_filt <- samp[samp$ffFam %in% match_df$ffFam,]

samp_9_filt <- samp_filt %>% filter(Group == "FC09")
samp_15_filt <- samp_filt %>% filter(Group == "FC15")


beta2 <- beta[,colnames(beta) %in% samp_filt$MethID]
beta_9 <- beta[,colnames(beta) %in% samp_9_filt$MethID]
beta_15 <- beta[,colnames(beta) %in% samp_15_filt$MethID]

# converting from betas to M-values
m <- logit(beta2)
m_9 <- logit(beta_9)
m_15 <- logit(beta_15)

# keeping only samples in INC file that were matched
INC2 <- INC[,colnames(INC) %in% colnames(m)]
INC2df <- as.data.frame(INC2)
INC2df <- INC2df[,match(colnames(m),colnames(INC2))]
INC2 <- as.matrix(INC2df)

# RUVM

# I originally had compared controls and exposed at age 15 and controls and exposed at age 9
# this is not what I used in my final analysis
# but I have left the code in case needed for reference ever
# comparison I actually used in my EWAS is below (not commented)
# compared year 9 to year 15 data in controls and repeated separately in exposed individuals


# from tutorial: "add negative control data to m-values"
#m_controls <- rbind(m, INC2)

# from tutorial: "create vector marking negative controls in data matrix"
#controls_vector <- rownames(m_controls) %in% rownames(INC2)

# making year 9 INC
# making sure samples are in same order as in m_9 file
# INC_9 <- INC2[,colnames(INC2) %in% colnames(m_9)]
# INC_9df <- as.data.frame(INC_9)
# INC_9df <- INC_9df[,match(colnames(m_9), colnames(INC_9df))]
# INC_9 <- as.matrix(INC_9df)

# creating one matrix with m values for samples and ther INCs at yr 9
# m_controls_9 <- rbind(m_9, INC_9)
# m_controls_9 <- rbind(m_9, INC_9)
#making vector marking with rows are actual cpg sites and which are controls
# controls_vector_9 <- rownames(m_controls_9) %in% rownames(INC_9)

# making year 15 INC
# # making sure samples are in same order as in m_15 file
# INC_15 <- INC2[,colnames(INC2) %in% colnames(m_15)]
# INC_15df <- as.data.frame(INC_15)
# INC_15df <- INC_15df[,match(colnames(m_15), colnames(INC_15df))]
# INC_15 <- as.matrix(INC_15df)

# m_controls_15 <- rbind(m_15, INC_15)
#making vector marking with rows are actual cpg sites and which are controls
# controls_vector_15 <- rownames(m_controls_15) %in% rownames(INC_15)

#exposure_vect <- c()
# making df of methid to exposure value
#for (i in samp_filt$MethID2){
#samp_filt2 <- samp_filt %>% filter(MethID2 == i)
#ffFam_1 <- as.character(samp_filt2$ffFam)
#match_df_filt <- match_df %>% filter(ffFam == ffFam_1)
#exp <- as.integer(match_df_filt$exposure) #need to change this if exp isn't binary
#exposure_vect <- c(exposure_vect, exp)
#}


# making dataframe linking methIDs to exposure value (0 or 1) 
# making sure this dataframe's methIDs in order of m value matrix
# exposure_vect_9 <- c()
# exposure_vect_9 <- c()
# for (i in colnames(m_controls_9)){
# samp_filt2 <- samp_filt %>% filter(MethID == i)
# ffFam_1 <- as.character(samp_filt2$ffFam)
# match_df_filt <- match_df %>% filter(ffFam == ffFam_1)
# exp <- as.integer(match_df_filt$exposure)
# exposure_vect_9 <- c(exposure_vect_9, exp)
# }


#making dataframe linking methIDs to exposure value (0 or 1) 
# making sure this dataframe's methIDs in order of m value matrix
# exposure_vect_15 <- c()
# for (p in colnames(m_controls_15)){
# samp_filt2 <- samp_filt %>% filter(MethID == p)
# ffFam_1 <- as.character(samp_filt2$ffFam)
# match_df_filt <- match_df %>% filter(ffFam == ffFam_1)
# exp <- as.integer(match_df_filt$exposure)
# exposure_vect_15 <- c(exposure_vect_15, exp)
# }



# df_exp_9 <- data.frame(colnames(m_controls_9), exposure_vect_9)
# df_exp_15 <- data.frame(colnames(m_controls_15), exposure_vect_15)
# 
# colnames(df_exp_9) <- c("MethID", "exposure")
# 
# colnames(df_exp_15) <- c("MethID", "exposure")

# vector that says whether a sample is exposed or not
# fact_exp_9 <- factor(df_exp_9$exposure, labels=c(0,1))
# fact_exp_15 <- factor(df_exp_15$exposure, labels=c(0,1))

#RUVM for age 9
# rfit1_9 <- RUVfit(Y=m_controls_9, X=fact_exp_9, ctl=controls_vector_9)
# rfit2_9 <- RUVadj(Y=m_controls_9, fit=rfit1_9)
# top1_9 <- topRUV(rfit2_9, num=Inf, p.BH = 1)
# 
# # next line finds CpGs least associated with factor of interest?
# ctl2_9 <- rownames(m_9) %in% rownames(top1_9[top1_9$p.BH_X1.1 > 0.5,])
# rfit3_9 <- RUVfit(Y=m_9, X=fact_exp_9, ctl=ctl2_9)
# rfit4_9 <- RUVadj(Y=m_9, fit=rfit3_9)
# 
# dmp_9 <- topRUV(rfit4_9)
# 
# cpg_9 <- c(cpg_9, rownames(dmp_9))
# F.p_9 <- c(F.p_9, dmp_9$F.p)
# F.p.BH_9 <- c(F.p.BH_9, dmp_9$F.p.BH)
# p_X1.1_9 <- c(p_X1.1_9, dmp_9$p_X1.1)
# p.BH_X1.1_9 <- c(p.BH_X1.1_9, dmp_9$p.BH_X1.1)
# b_X1.1_9 <- c(b_X1.1_9, dmp_9$b_X1.1)
# sigma2_9 <- c(sigma2_9, dmp_9$sigma2)
# var.b_X1.1_9 <- c(var.b_X1.1_9, dmp_9$var.b_X1.1)
# fit.ctl_9 <- c(fit.ctl_9, dmp_9$fit.ctl)
# mean_9 <- c(mean_9, dmp_9$mean)

#RUVM for age 15
# rfit1_15 <- RUVfit(Y=m_controls_15, X=fact_exp_15, ctl=controls_vector_15)
# rfit2_15 <- RUVadj(Y=m_controls_15, fit=rfit1_15)
# top1_15 <- topRUV(rfit2_15, num=Inf, p.BH=1)
# ctl2_15 <- rownames(m_15) %in% rownames(top1_15[top1_15$p.BH_X1.1 > 0.5,])
# rfit3_15 <- RUVfit(Y=m_15, X=fact_exp_15, ctl=ctl2_15)
# rfit4_15 <- RUVadj(Y=m_15, fit=rfit3_15)
# 
# dmp_15 <- topRUV(rfit4_15)
# 
# cpg <- c(cpg, rownames(dmp_15))
# F.p <- c(F.p, dmp_15$F.p)
# F.p.BH <- c(F.p.BH, dmp_15$F.p.BH)
# p_X1.1 <- c(p_X1.1, dmp_15$p_X1.1)
# p.BH_X1.1 <- c(p.BH_X1.1, dmp_15$p.BH_X1.1)
# b_X1.1 <- c(b_X1.1, dmp_15$b_X1.1)
# sigma2 <- c(sigma2, dmp_15$sigma2)
# var.b_X1.1 <- c(var.b_X1.1, dmp_15$var.b_X1.1)
# fit.ctl <- c(fit.ctl, dmp_15$fit.ctl)
# mean <- c(mean, dmp_15$mean)



#comparing year 9 and 15 in controls
#control comparison in loop but exposed comparison out of loop because
# only done once for exposed
# keeping only control samples
samp_c <- samp[samp$ffFam %in% control_ffFam,]
beta_c <- beta[,colnames(beta) %in% samp_c$MethID]
#converting to m values (from betas)
m_c <- logit(beta_c)
#making INC file that only includes matched control samples
INC2_c <- INC[,colnames(INC) %in% colnames(m_c)]
# making sure file has samples in same order as M file
INC2df_c <- as.data.frame(INC2_c)
INC2df_c <- INC2df_c[,match(colnames(m_c),colnames(INC2_c))]
INC2_c <- as.matrix(INC2df_c)

# making combined df with INC and m values
m_controls_c <- rbind(m_c, INC2_c)

# vector to denote whether a methID is year 15 or year 9
year_vect_c <- c()

for (i in colnames(m_controls_c)){
samp_filt2c <- samp_c %>% filter(MethID == i)
year1 <- samp_filt2c$Group
year_vect_c <- c(year_vect_c, year1)
}

year_vect_c <- ifelse(year_vect_c=="FC09",0,1)
fact_year_vect_c <- factor(year_vect_c, labels=c(0,1))

# marking which rows are CpG sites and which are controls
controls_vector_c <- rownames(m_controls_c) %in% rownames(INC2_c)

#RUVM
rfit1_c <- RUVfit(Y=m_controls_c, X=fact_year_vect_c, ctl=controls_vector_c)
rfit2_c <- RUVadj(Y=m_controls_c, fit=rfit1_c)
top1_c <- topRUV(rfit2_c, num=Inf)

# next line finds CpGs least associated with factor of interest?
ctl2_c <- rownames(m_c) %in% rownames(top1_c[top1_c$p.BH_X1.1 > 0.5,])
rfit3_c <- RUVfit(Y=m_c, X=fact_year_vect_c, ctl=ctl2_c)
rfit4_c <- RUVadj(Y=m_c, fit=rfit3_c)

dmp_c <- topRUV(rfit4_c, num=Inf)



cpg_c <- c(cpg_c, rownames(dmp_c))
F.p_c <- c(F.p_c, dmp_c$F.p)
F.p.BH_c <- c(F.p.BH_c, dmp_c$F.p.BH)
p_X1.1_c <- c(p_X1.1_c, dmp_c$p_X1.1)
p.BH_X1.1_c <- c(p.BH_X1.1_c, dmp_c$p.BH_X1.1)
b_X1.1_c <- c(b_X1.1_c, dmp_c$b_X1.1)
sigma2_c <- c(sigma2_c, dmp_c$sigma2)
var.b_X1.1_c <- c(var.b_X1.1_c, dmp_c$var.b_X1.1)
fit.ctl_c <- c(fit.ctl_c, dmp_c$fit.ctl)
mean_c <- c(mean_c, dmp_c$mean)


}

# df <- data.frame(cpg,F.p, F.p.BH, p_X1.1, p.BH_X1.1, b_X1.1, sigma2, var.b_X1.1, fit.ctl, mean) 
# df_9 <- data.frame(cpg_9, F.p_9, F.p.BH_9, p_X1.1_9, p.BH_X1.1_9, b_X1.1_9, sigma2_9, var.b_X1.1_9, fit.ctl_9, mean_9)

# control df
df_c <- data.frame(cpg_c, F.p_c, F.p.BH_c, p_X1.1_c, p.BH_X1.1_c, b_X1.1_c, sigma2_c, var.b_X1.1_c, fit.ctl_c, mean_c)



# repeating procedure comparing year 9 and year 15 data to find dif methylated sites 
# except this time doing this for exposed group (just 1 iteration because just 1 exposed group)

# keeping just exposed individuals
samp_exp <- samp[samp$ffFam %in% exposed_ffFam,]
beta_exp <- beta[,colnames(beta) %in% samp_exp$MethID]
m_exp <- logit(beta_exp)
#making INC file that only includes matched control samples
INC2_exp <- INC[,colnames(INC) %in% colnames(m_exp)]
INC2df_exp <- as.data.frame(INC2_exp)
# making sure file has samples in same order as M file
INC2df_exp <- INC2df_exp[,match(colnames(m_exp),colnames(INC2_exp))]
INC2_exp <- as.matrix(INC2df_exp)

# making combined df with INC and m values
m_controls_exp <- rbind(m_exp, INC2_exp)

# vector to denote whether a methID is year 15 or year 9
year_vect <- c()

for (i in colnames(m_controls_exp)){
samp_filt2 <- samp_exp %>% filter(MethID == i)
year1 <- samp_filt2$Group
year_vect <- c(year_vect, year1)
}

year_vect <- ifelse(year_vect=="FC09",0,1)
fact_year_vect <- factor(year_vect, labels=c(0,1))

# marking which rows are CpG sites and which are controls
controls_vector_exp <- rownames(m_controls_exp) %in% rownames(INC2_exp)

#RUVM
rfit1_exp <- RUVfit(Y=m_controls_exp, X=fact_year_vect, ctl=controls_vector_exp)
rfit2_exp <- RUVadj(Y=m_controls_exp, fit=rfit1_exp)
top1_exp <- topRUV(rfit2_exp, num=Inf)

# next line finds CpGs least associated with factor of interest?
ctl2_exp <- rownames(m_exp) %in% rownames(top1_exp[top1_exp$p.BH_X1.1 > 0.5,])
rfit3_exp <- RUVfit(Y=m_exp, X=fact_year_vect, ctl=ctl2_exp)
rfit4_exp <- RUVadj(Y=m_exp, fit=rfit3_exp)

dmp_exp <- topRUV(rfit4_exp, number=Inf)


# these lines are left in from when I had originally compared exposed to controls at age 15 and at age 9
# instead of, now, comparing controls age 9 and 15 and separately comparing exposed at ages 9 and 15
#write.csv(df, "/home/HPA/Hannah/k450_meth_2-14_100_15_fixed.csv")
#write.csv(df_9, "/home/HPA/Hannah/k450_meth_2-14_100_9_fixed.csv")


write.csv(dmp_exp, "/home/HPA/Hannah/k450_meth_2-14_exp_year915.csv")
write.csv(df_c, "/home/HPA/Hannah/k450_meth_2-14_c_year915.csv")
