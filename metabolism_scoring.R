#metabolism scoring
#code name on server: /home/HPA/HPA_group_datasets/ADME/Hannah/scoring_updated_3-5.R

library(dplyr)
library(tidyr)
library(stringr)
library(vcfR)

#scoring all ffFam IDs in dist

# reading in fracking distance file
dist <- readRDS("/home/HPA/Hannah/BenCode/datafiles/frackDist.RDS")

# reading in file of metabolism variants
var <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/metabolismvariants_2-21_4.csv")

# reading in VCF file
vcf <- read.vcfR("/home/HPA/HPA_group_datasets/ADME/Hannah/hannah_var_regions_extracted/topmed.adme.vcf.gz")

# using VCFR package to put VCF into tidyr object format
vcf1 <- vcfR2tidy(vcf, format_fields=c("GT", "DP"))
vcf1 <- vcf1$gt

#fixing chr
chr_fix <- c()
for (i in vcf1$ChromKey){
if (i==1){chr_fix <- c(chr_fix, 10)}
if (i==2){chr_fix <- c(chr_fix, 15)}
if (i==3){chr_fix <- c(chr_fix, 19)}
if (i==4){chr_fix <- c(chr_fix, 4)}
if (i==5){chr_fix <- c(chr_fix, 7)}
}
vcf1$ChromKey <- chr_fix

#splitting alleles into two columns
allele1 <- str_sub(vcf1$gt_GT_alleles, end=-3)
allele2 <- str_sub(vcf1$gt_GT_alleles, start=3)
vcf1$allele1 <- allele1
vcf1$allele2 <- allele2

#removing any NAs
vcf1 <- vcf1[!is.na(vcf1$gt_GT),]
vcf1 <- vcf1[!is.na(vcf1$gt_GT_alleles),]

#making combined position column
vcf1$combined <- paste(vcf1$ChromKey, vcf1$POS, sep=":")

# regions=regions of variants from variants file to use (other variants weren't available in any/most individuals)
regions <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/varpos_numindiv")
colnames(regions) <- c("row", "pos", "counts")
regions_pos <- regions$pos[regions$counts >= 2800]

vcf1 <- vcf1[vcf1$combined %in% regions_pos,]
ffid <- read.csv("/home/HPA/HPA_group_datasets/ADME/Hannah/ffid_19var")

#filtering individuals who have all of the 12 variants
vcf1 <- vcf1[vcf1$Indiv %in% ffid$x,]



#for single loci decreasing variants
vcf1_1c <- vcf1 %>% filter(combined == "15:74745879")
vcf1_1c$score <- str_count(vcf1_1c$gt_GT_alleles, "A")

vcf1_2e16 <- vcf1 %>% filter(combined == "10:133535040")
vcf1_2e16$score <- str_count(vcf1_2e16$gt_GT_alleles, "A")

vcf1_2a618 <- vcf1 %>% filter(combined == "19:40844759")
vcf1_2a618$score <- str_count(vcf1_2a618$gt_GT_alleles, "T")

vcf1_3a4 <- vcf1 %>% filter(combined == "7:99768693")
vcf1_3a4$score <- str_count(vcf1_3a4$gt_GT_alleles, "T")

vcf1_adh1c <- vcf1 %>% filter(combined == "4:99342808")
vcf1_adh1c$score <- str_count(vcf1_adh1c$gt_GT_alleles, "T")

#for single loci increasing variants
vcf1_adh1b <- vcf1 %>% filter(combined == "4:99307860")
vcf1_adh1b$score <- (str_count(vcf1_adh1b$gt_GT_alleles, "T") * -1)

# for 2 loci decreasing variants
vcf2e15 <- vcf1 %>% filter(combined == "10:133526101" | combined=="10:133526341")



#making reshaped df (one row per individual with columns for each locus of variant)
e15_a1_pos1 <- c()
e15_a2_pos1 <- c()
e15_a1_pos2 <- c()
e15_a2_pos2 <- c()

vcf2e15_pos1 <- vcf2e15 %>% filter(combined == "10:133526101")
vcf2e15_pos2 <- vcf2e15 %>% filter(combined == "10:133526341")
for (j in unique(vcf2e15$Indiv)){
vcf2e15_pos1filt <- vcf2e15_pos1 %>% filter(Indiv == j)
vcf2e15_pos2filt <- vcf2e15_pos2 %>% filter(Indiv == j)
e15_a1_pos1 <- c(e15_a1_pos1, vcf2e15_pos1filt$allele1[1])
e15_a2_pos1 <- c(e15_a2_pos1, vcf2e15_pos1filt$allele2[1])
e15_a1_pos2 <- c(e15_a1_pos2, vcf2e15_pos2filt$allele1[1])
e15_a2_pos2 <- c(e15_a2_pos2, vcf2e15_pos2filt$allele2[1])
}
e15_df <- data.frame(unique(vcf2e15$Indiv), e15_a1_pos1, e15_a2_pos1, e15_a1_pos2, e15_a2_pos2)
colnames(e15_df) <- c("Indiv", "a1pos1", "a2pos1", "a1pos2", "a2pos2")
e15_df$score1 <- ifelse(e15_df$a1pos1 == "C" & e15_df$a1pos2 == "T", 1, 0)
e15_df$score2 <- ifelse(e15_df$a2pos1 == "C" & e15_df$a2pos2 == "T", 1, 0)
e15_df$score <- e15_df$score1 + e15_df$score2




# for 2b6 variants (different because one variant is LOF)

#for single loci increasing variants
vcf1_2b622 <- vcf1 %>% filter(combined == "19:40991224")
vcf1_2b622$score <- (str_count(vcf1_2b622$gt_GT_alleles, "C") * -1)

#for lof variant
#for single loci increasing variants
vcf1_2b618 <- vcf1 %>% filter(combined == "19:41012316")
vcf1_2b618$score <- str_count(vcf1_2b618$gt_GT_alleles, "C")

score_2b6 <- c()
for (p in unique(vcf1_2b622$Indiv)){
vcf1_2b618_filt <- vcf1_2b618 %>% filter(Indiv == p)
vcf1_2b622_filt <- vcf1_2b622 %>% filter(Indiv == p)
score_18 <- vcf1_2b618$score[1] 
score_22 <- vcf1_2b622$score[1]

score_tot <- score_18 + score_22
if(vcf1_2b618_filt$score[1] == 2){
score_2b6 <- c(score_2b6, vcf1_2b618_filt$score[1])
}
else {score_2b6 <- c(score_2b6, score_tot)}
}
df_2b6 <- data.frame(unique(vcf1_2b622$Indiv), score_2b6)
colnames(df_2b6) <- c("Indiv", "score")


#extracting scores for each variant
vcf1_1c <- vcf1_1c %>% select(Indiv, score)
vcf1_2e16 <- vcf1_2e16 %>% select(Indiv, score)
vcf1_2a618 <- vcf1_2a618 %>% select(Indiv, score)
vcf1_3a4 <- vcf1_3a4 %>% select(Indiv, score)
vcf1_adh1c <- vcf1_adh1c %>% select(Indiv, score)
vcf1_adh1b <- vcf1_adh1b %>% select(Indiv, score)
e15_df <- e15_df %>% select(Indiv, score)
k1_df <- k1_df %>% select(Indiv, score)
vcf1_2a623 <- vcf1_2a623 %>% select(Indiv, score)
df_2b6 <- df_2b6 %>% select(Indiv, score)

# keeping only individuals with scores for each variant
intersect_indiv <- Reduce(intersect, list(vcf1_1c$Indiv, vcf1_2e16$Indiv, vcf1_2a618$Indiv,vcf1_3a4$Indiv, vcf1_adh1c$Indiv, vcf1_adh1b$Indiv, df_2b6$Indiv, e15_df$Indiv))

#renaming columns before merging and filtering to only have individuals that have scores for all 9 var
vcf1_1c <- vcf1_1c[vcf1_1c$Indiv %in% intersect_indiv,]
colnames(vcf1_1c) <- c("Indiv", "score_1c")
vcf1_2e16 <- vcf1_2e16[vcf1_2e16$Indiv %in% intersect_indiv,]
colnames(vcf1_2e16) <- c("Indiv", "score_2e16")
vcf1_2a618 <- vcf1_2a618[vcf1_2a618$Indiv %in% intersect_indiv,]
colnames(vcf1_2a618) <- c("Indiv", "score_2a618")
vcf1_3a4 <- vcf1_3a4[vcf1_3a4$Indiv %in% intersect_indiv,]
colnames(vcf1_3a4) <- c("Indiv", "score_3a4")
vcf1_adh1c <- vcf1_adh1c[vcf1_adh1c$Indiv %in% intersect_indiv,]
colnames(vcf1_adh1c) <- c("Indiv", "score_adh1c")
vcf1_adh1b <- vcf1_adh1b[vcf1_adh1b$Indiv %in% intersect_indiv,]
colnames(vcf1_adh1b) <- c("Indiv", "score_adh1b")
e15_df <- e15_df[e15_df$Indiv %in% intersect_indiv,]
colnames(e15_df) <- c("Indiv", "score_e15")
df_2b6 <- df_2b6[df_2b6$Indiv %in% intersect_indiv,]
colnames(df_2b6) <- c("Indiv", "score_2b6")

#merging to make one df with all scores for an individual
df_scores <- merge(vcf1_1c, vcf1_2e16, by="Indiv")
df_scores <- merge(df_scores, vcf1_2a618, by="Indiv")
df_scores <- merge(df_scores, vcf1_3a4, by="Indiv")
df_scores <- merge(df_scores, vcf1_adh1c, by="Indiv")
df_scores <- merge(df_scores, vcf1_adh1b, by="Indiv")
df_scores <- merge(df_scores, e15_df, by="Indiv")
df_scores <- merge(df_scores, df_2b6, by="Indiv")

# creating total score for a person
df_scores <- df_scores %>% group_by(Indiv) %>% mutate(total = sum(score_1c, score_2e16, score_2a618, score_3a4, score_adh1c, score_adh1b, score_e15, score_2b6))
