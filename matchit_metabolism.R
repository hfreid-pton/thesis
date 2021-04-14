## new matchit code for matching during metabolism analysis
# code name on server: /home/HPA/Hannah/matchit_11-1_metabolism.R


library(dplyr)
library(MatchIt)
library(tidyr)
library(stringr)


#read in score data and create binary scores
score <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/9var_scores_3-14.csv")
score$ffFam <- str_sub(score$Indiv, end=-5)
score$binary <- ifelse(score$total < 5,0,1)
score$total_alc <- (score$score_2e16 + score$score_adh1c + score$ score_adh1b + score$score_e15)
score$binary_alc <- ifelse(score$total_alc < 2,0,1)



##these lines with data and processing for epigenetic clock analysis are uncommented when 
# using this code to match for epigenetic clock metabolism analyses
# but stay commented when using this code for health metablism analyses
#read in epigenetic clock data
#clocks <- read.csv("/home/HPA/HPA_group_datasets/MethylationData/January2021/epic_epigeneticAges.csv")

#read in problem sample ids to remove problem samples before matching
#prob1 <- readRDS("/home/HPA/Hannah/BenCode/datafiles/probIDsEpic.RDS")
#prob2 <- read.csv("/home/HPA/Hannah/prob_allEPIC_hannah.csv", stringsAsFactors=FALSE)
#prob3 <- read.csv("/home/HPA/Hannah/prob_clock_allepic_hannah.csv", stringsAsFactors=FALSE)
#prob <- c(prob1, prob2$prob_both, prob3$x)
#prob <- unique(prob)

#ffFam <- str_sub(clocks$ffID, end=-5)
#clocks$ffFam <- ffFam


#removing problem samples
#clocks <- clocks[!clocks$ffFam %in% prob,]


# keeping only people with scores and clock data
#clocks <- clocks[clocks$ffFam %in% score$ffFam,]


#keeping only people with slow metabolism scores
# this line is commented out when doing the comparison within the exposed group
# but used when using this for matching slow exp to slow controls
score_slow <- score$ffFam[score$binary == 1] #this line is changed to score$binary_alc==1 for alc score
#clocks <- clocks[clocks$ffFam %in% score_slow,]



# reading in distance from fracking site data
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")

# defining exposed group
dist_10 <- dist %>% filter(distance_ff_site <= 10)
all_ffid <- unique(dist$ffid)


ffid_10 <- dist_10$ffid
counts_10 <- c()

for (i in 1:length(all_ffid)){
num_10 <- length(ffid_10[ffid_10 == all_ffid[i]])
counts_10 <- append(counts_10, num_10)
}

num_wells <- data.frame(all_ffid, counts_10)


library(stringr)
ffam_test <- str_sub(num_wells$all_ffid, end=-5)
num_wells$ffFam <- ffam_test

# reading in survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

exposed <- num_wells %>% filter(counts_10 > 0)


# defining potential controls
within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam


library(stringr)


num_wells_10 <- num_wells
num_wells_10$exposure <- ifelse(num_wells_10$counts_10 > 0, 1, 0)
num_wells_10 <- num_wells_10[(num_wells_10$ffFam %in% controls_ffFam) | (num_wells_10$ffFam %in% exposed$ffFam),]

merged2 <- surv[surv$ffFam %in% num_wells_10$ffFam,]

merged3 <- merge(merged2, num_wells_10, by="ffFam")

# dropping individuals with NAs or error codes for covariates used in matching
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

#only keeping individuals who either didn't move or who moved but didn't switch schools
merged4 <- merged4 %>% filter(p6j1 == 2 | p6j5 == 0)


merged4 <- merged4 %>% select(ffFam, cm1ethrace,exposure,cm1bsex, cm5povca, p5h15c, m1g4,k6d46, ck6pcgrel, cp6povca)



#making cm1ethrace into dummary var
merged4$his <- ifelse(merged4$cm1ethrace == 3, 1, 0)
merged4$afr <- ifelse(merged4$cm1ethrace == 2, 1, 0)


# this line is uncommented when this code is used in comparisons of slow exp to fast exp
# this line is here to just keep individuals with scores in
#merged4 <- merged4[merged4$ffFam %in% score$ffFam,]


# matching
m.out <- matchit(exposure ~ his + afr + cm1bsex + cm5povca + p5h15c + m1g4 + k6d46 + ck6pcgrel + cp6povca, data=merged4, method="nearest", ratio=1)







