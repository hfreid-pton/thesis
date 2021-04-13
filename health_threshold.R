# finding distance threshold of fracking effects
# server code name: /home/HPA/Hannah/health_binary_difdist_rem_5_2.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(MatchIt)
library(stringr)

dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")

# radius is upper limit of range 
radius <- seq(from=0, to=60, by =5)
p_p6b19 <- c() #seeing where health effects (p6b19 and p6b23) are no longer assoc. w fracking exposure
p_p6b23 <- c()
num_exposed <- c()
upper_lim <- c()
ffFam_exposed_used <- c()
#number of exposed people just at that radius (other num_exposed is for all individuals in closer radii too)
num_exposed2 <- c() 


# running loop--go through loop once for each radius range (i.e 0-5 km, 5-10 km)
# each radius-range group of exposed individuals matched to one set of controls
# then compare to see association between fracking exposure (at certain distance) and 2 health cond.
for (i in 1:length(radius)){
upper_lim_1 = radius[(i+1)]
upper_lim <- c(upper_lim, upper_lim_1)
dist_rad <- dist %>% filter(distance_ff_site < upper_lim_1 & distance_ff_site >= radius[i])
all_ffid <- unique(dist$ffid)

ffid_rad <- dist_rad$ffid

counts_rad <- c()

# getting number of wells each ffid is near within the radius
for (n in 1:length(all_ffid)){
num_rad <- length(ffid_rad[ffid_rad == all_ffid[n]])
counts_rad <- append(counts_rad, num_rad)
}

num_wells <- data.frame(all_ffid, counts_rad)

# converting ffid to ffFam
ffam_test <- str_sub(num_wells$all_ffid, end=-5)
num_wells$ffFam <- ffam_test


#reading in survey data
surv <- read.csv("/home/HPA/HPA_group_datasets/ffsurvey/FF_allwaves_2018_FFID.csv")

#defining exposed individuals at each radius vs. controls (always more than 60 km away)
exposed <- num_wells %>% filter(counts_rad > 0)


#finding potential controls to be used for matching
within_60 <- dist %>% filter(distance_ff_site <= 60)
controls <- setdiff(dist$ffid, within_60$ffid)
controls_ffFam <- str_sub(controls, end=-5)

exposed_ffFam <- exposed$ffFam



num_wells_rad <- num_wells
num_wells_rad$exposure <- ifelse(num_wells_rad$counts_rad > 0, 1, 0)
num_wells_rad <- num_wells_rad[(num_wells_rad$ffFam %in% controls_ffFam) | (num_wells_rad$ffFam %in% exposed$ffFam),]

# removing ffFam already used so those exposed are only exposed at that range
num_wells_rad <- num_wells_rad %>% filter(!ffFam %in% ffFam_exposed_used)
num_wells_rad_exposed_temp <- num_wells_rad %>% filter(exposure == 1)
ffFam_exposed_used <- c(ffFam_exposed_used, num_wells_rad_exposed_temp$ffFam)
num_exposed <- c(num_exposed, length(num_wells_rad_exposed_temp$ffFam))


merged2 <- surv[surv$ffFam %in% num_wells_rad$ffFam,]

merged3 <- merge(merged2, num_wells_rad, by="ffFam")

#dropping NAs/error codes for variables used in matching
merged4 <- drop_na(merged3,c("cm1ethrace", "cm1bsex", "cm5povca", "p5h15c", "m1g4", "k6d46", "ck6pcgrel", "cp6povca"))

merged4 <- merged4[merged4$cm1bsex > 0,]
merged4 <- merged4[merged4$cm5povca > 0,]
merged4 <- merged4[merged4$p5h15c >= 0,]
merged4 <- merged4[merged4$m1g4 > 0,]
merged4 <- merged4[merged4$k6d46 >= 0,]
merged4 <- merged4[merged4$ck6pcgrel > 0,]
merged4 <- merged4[merged4$cp6povca > 0,]
merged4 <- merged4[merged4$cm1ethrace > 0,]

#keeping only people who didn't move or who moved and didn't switch schools
merged4 <- merged4 %>% filter(p6j1 == 2 | p6j5 == 0)

merged4 <- merged4 %>% select(ffFam, cm1ethrace, other, exposure,cm1bsex, cm5povca, p5h15c, m1g4,k6d46, ck6pcgrel, cp6povca)

# making cm1ethrace a dummy var (white=reference)
merged4$his <- ifelse(merged4$cm1ethrace == 3, 1, 0)
merged4$afr <- ifelse(merged4$cm1ethrace == 2, 1, 0)

# matching to one iteration of controls
m.out <- matchit(exposure ~ his + afr + other + cm1bsex + cm5povca + p5h15c + m1g4 + k6d46 + ck6pcgrel + cp6povca, data=merged4, method="nearest", ratio=1)





match_df <- match.data(m.out)

num_exposed2 <- c(num_exposed2, length(intersect(num_wells_rad_exposed_temp$ffFam, match_df$ffFam)))


num_wells_rad <- match_df



#trouble breathing
p6b19 <- c()

# seeing dr for illness in past year
p6b23 <- c()



for (j in 1:length(num_wells_rad$ffFam)){
merged_i <- merged[merged$ffFam == num_wells_rad$ffFam[j], ]
p6b19_1 <- as.integer(merged_i$p6b19[1])
p6b23_1 <- as.integer(merged_i$p6b23[1])
p6b19 <- append(p6b19, p6b19_1)
p6b23 <- append(p6b23, p6b23_1)
}



# p6b19 (same as health_analysis.R)
num_wells_rad$p6b19 <- p6b19
num_wells_rad_p6b19 <- num_wells_rad %>% filter(p6b19 >= 1)
num_wells_rad_p6b19$p6b19_corrected <- ifelse(num_wells_rad_p6b19$p6b19 == 2, 0, 1)

num_wells_rad_p6b19_expo <- num_wells_rad_p6b19 %>% filter(exposure == 1)
num_wells_rad_p6b19_c <- num_wells_rad_p6b19 %>% filter(exposure == 0)
num_wells_rad_p6b19_expo <- num_wells_rad_p6b19_expo$p6b19_corrected[is.na(num_wells_rad_p6b19_expo$p6b19_corrected) == F]
num_wells_rad_p6b19_c <- num_wells_rad_p6b19_c$p6b19_corrected[is.na(num_wells_rad_p6b19_c$p6b19_corrected) == F]
exp_p6b19_rad <- c(length(num_wells_rad_p6b19_expo[num_wells_rad_p6b19_expo == 1]), length(num_wells_rad_p6b19_expo[num_wells_rad_p6b19_expo == 0]))
c_p6b19_rad <- c(length(num_wells_rad_p6b19_c[num_wells_rad_p6b19_c == 1]), length(num_wells_rad_p6b19_c[num_wells_rad_p6b19_c == 0]))

counts_p6b19_corrected <- chisq.test(as.table(rbind(exp_p6b19_rad, c_p6b19_rad)))

counts_p6b19_notcorrected <- chisq.test(as.table(rbind(exp_p6b19_rad, c_p6b19_rad)), correct=F)


p6b19_z <- prop.test(x=c(exp_p6b19_rad[1], c_p6b19_rad[1]), n=c(sum(exp_p6b19_rad), sum(c_p6b19_rad)))

p_p6b19 <- c(p_p6b19, prop_p6b19_notcorrected$p.value)


#p6b23
num_wells_rad$p6b23 <- p6b23
num_wells_rad_p6b23 <- num_wells_rad %>% filter(p6b23 >= 1)
num_wells_rad_p6b23$p6b23_corrected <- ifelse(num_wells_rad_p6b23$p6b23 == 2, 0, 1)
num_wells_rad_p6b23_expo <- num_wells_rad_p6b23 %>% filter(exposure == 1)
num_wells_rad_p6b23_c <- num_wells_rad_p6b23 %>% filter(exposure == 0)
num_wells_rad_p6b23_expo <- num_wells_rad_p6b23_expo$p6b23_corrected[is.na(num_wells_rad_p6b23_expo$p6b23_corrected) == F]
num_wells_rad_p6b23_c <- num_wells_rad_p6b23_c$p6b23_corrected[is.na(num_wells_rad_p6b23_c$p6b23_corrected) == F]
exp_p6b23_rad <- c(length(num_wells_rad_p6b23_expo[num_wells_rad_p6b23_expo == 1]), length(num_wells_rad_p6b23_expo[num_wells_rad_p6b23_expo == 0]))
c_p6b23_rad <- c(length(num_wells_rad_p6b23_c[num_wells_rad_p6b23_c == 1]), length(num_wells_rad_p6b23_c[num_wells_rad_p6b23_c == 0]))

counts_p6b23_corrected <- chisq.test(as.table(rbind(exp_p6b23_rad, c_p6b23_rad)))
counts_p6b23_notcorrected <- chisq.test(as.table(rbind(exp_p6b23_rad, c_p6b23_rad)), correct=F)


p_p6b23 <- c(p_p6b23, prop_p6b23_notcorrected$p.value)



print(i)

}


df <- data.frame(lower_lim, radius, num_exposed, p_p6b19, p_p6b23)

plot_p6b19 <- ggplot(data=df, mapping=aes(x=radius, y=p_p6b19))
plot_p6b19 + geom_point()




p <- ggplot()
p <- p + geom_line(data=df, aes(x=radius, y=p_p6b19, color="Chest/Breathing Problems"))
p <- p + geom_line(data=df, aes(x=radius, y=p_p6b23, color="Seeing Doctor for Illness"))
p <- p + labs(x="Radius", y="P-value for fracking exposure") + scale_x_continuous(breaks=seq(0,40,5)) + scale_y_continuous(breaks=seq(0,1,0.1))

