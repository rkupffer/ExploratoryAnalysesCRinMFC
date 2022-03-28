###################################################  Exploratory Analyses  ###################################################
load("5_rq2_lpa_results.RData")
all_items <- read.csv2("2_All_Items_Coding.csv", sep = ",")[,2:8]
df <- readRDS("4_indices.rds")
bin <- readRDS("2_df_recoded.rds")
bin <- subset(bin, select = c(ID, b.bi.1:b.bi.60))

library(CRinMFC)
library(multiplex)
library(psych)
library(apaTables)
library(tidyLPA)
library(MFCblockInfo)


est_lpa <- get_estimates(lpa)
est_mod2.4 <- as.data.frame(lpa$model_2_class_4$dff)


# additional functions which are not part of any package
source("0_add_functions.R")

#_________________________________ age and gender variables ______________________________ ####
# all demographic information was removed from due to data privacy 

names(df)[names(df) == "DA02_01"] <- "age"
names(df)[names(df) == "DA03"] <- "gender"
df$gender <- car::recode(df$gender, recodes = "1=0; 2=1; 3=NA; NA=NA")
df_0 <- data.frame(ID = df$ID, age = df$age, gender = df$gender)


#_________________________________ Criteria for being careless ______________________________ ####
# A: participants who were assigned to LPA-Class 2-4 (Model 3)
df$lpa_class <- est_mod2.4$Class
df$lpa_class <- car::recode(df$lpa_class, recodes = "1=2; 2=1; 3=4; 4=3") # but like this its in line with the order above
df$critA <- ifelse(df$lpa_class > 1, 1, 0)
prop.table(table(df$critA)) #check

# B: participants who met the criteria of being careless on the indices Triplet Variance, RTI, LCM, and Instructed Response Triplet
df$critB <- ifelse(df$rti == 1 | df$lcm == 1 | df$tv == 1 | df$ir != 0, 1 , 0)
df[is.na(df$critB), which(colnames(df)== "ir")] # 7 participants have a missing on one of the irt and have therefore also an missing on critB
df[, "critB"] <- car::recode(df[, "critB"], recodes = "NA=0")
prop.table(table(df$critB))

# A x B
table(df$critA, df$critB)
round2(prop.table(table(df$critA, df$critB)), 2)

#_________________________________ BIG FIVE INVENTORY 2 ______________________________ ####
#--------------------------------- prepare data ####
# ID, age, gender, bfi (60), critA, critB
datcriteria <- df[, c("ID", "critA", "critB")]
datvars <- df[ , c(which(colnames(df) == "ID"), 
                   which(colnames(df) == "age"), 
                   which(colnames(df) == "gender"))]
datTIRT_bi <- merge(datcriteria, datvars, by = "ID")
datTIRT_bi <- merge(datTIRT_bi, bin, by = "ID")
dim(datTIRT_bi)
head(datTIRT_bi, 3)

datTIRT_bi[,4:5] <- lapply(datTIRT_bi[,4:5], function(x) car::recode(x, recodes = "NA=-99"))
datTIRT_bi[,6:65] <- lapply(datTIRT_bi[,6:65], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))
datTIRT_bi$age <- datTIRT_bi$age/10
write.dat(datTIRT_bi, "6_exploratory")


#--------------------------------- critA ####
# prepare TIRT syntax
m_dload_bi <- designLoadings(all_items, quest = "BI", no.traits = 5)
colnames(m_dload_bi) <- c("INeu", "IExt", "IOpe", "IAgr", "ICon")
rownames(m_dload_bi) <- ifelse((1:nrow(m_dload_bi))<10,paste0("I","0",1:nrow(m_dload_bi)), paste0("I",1:nrow(m_dload_bi)))
tirt_bi_critA <- tirt.mplus.syntax.mod(design.load = m_dload_bi, 
                                       names.pairs = NULL,
                                       item.short = "I",
                                       id.var = "ID",
                                       file.data = "datTIRT_bi.dat", 
                                       title = "06_TIRT_BFI-critA", 
                                       out.command = "sampstat standardized;",
                                       fscores.file = "fscores_bi_critA.dat",
                                       missings.as = "-99", 
                                       all.var = c("critA", "critB", "age", "gender"),
                                       use.obs = c("critA EQ 0"))

# save as input-file
cat(paste(tirt_bi_critA, collapse="\n\n"), file="6_exploratory/06_TIRT_BFI-critA.inp")

#--------------------------------- critB ####
tirt_bi_critB <- tirt.mplus.syntax.mod(design.load = m_dload_bi, 
                                       names.pairs = NULL,
                                       item.short = "I",
                                       id.var = "ID",
                                       file.data = "datTIRT_bi.dat", 
                                       title = "06_TIRT_BFI-critB", 
                                       out.command = "sampstat standardized;",
                                       fscores.file = "fscores_bi_critB.dat",
                                       missings.as = "-99", 
                                       all.var = c("critA", "critB", "age", "gender"),
                                       use.obs = c("critB EQ 0"))

# save as input-file
cat(paste(tirt_bi_critB, collapse="\n\n"), file="6_exploratory/06_TIRT_BFI-critB.inp")


#_________________________________ Empirical Reliability ______________________________ ####

# compare the empirical reliability of trait estimates on the scales based on the whole sample 
# with that in the sample without the careless respondents

#--------------------------------- BFI #### 
#--------------------------------- all participants ####
m_theta.bfi_all <- read.table("3_Mplus/fscores_bi.dat", header = FALSE, na.strings = "*")[, c(71,61,63,65,67,69)]
colnames(m_theta.bfi_all)<- c("ID", paste0(c("N","E","O","A","C"), "_bi_all"))
m_ses.bi_all <- read.table("3_Mplus/fscores_bi.dat", header = FALSE, na.strings = "*")[,c(62,64,66,68,70)]
colnames(m_ses.bi_all)<- c(paste0(c("N","E","O","A","C"),"_se"))

df_bi <- merge(df_0, m_theta.bfi_all, by = "ID", all = TRUE)

(emp_rel_bfi_all <- MFCblockInfo::calc.rel.emp(m_theta.bfi_all[, 2:6], m_ses.bi_all))
(emp_rel_bfi_all_z <- psych::fisherz(sqrt(emp_rel_bfi_all)))

#--------------------------------- critA ####
m_theta.bfi_critA <- read.table("6_exploratory/fscores_bi_critA.dat", header = FALSE, na.strings = "*")[, c(71,61,63,65,67,69)]
colnames(m_theta.bfi_critA)<- c("ID", paste0(c("N","E","O","A","C"), "_bi_critA"))
m_ses.bi_critA  <- read.table("6_exploratory/fscores_bi_critA.dat", header = FALSE, na.strings = "*")[,c(62,64,66,68,70)]
colnames(m_ses.bi_critA )<- c(paste0(c("N","E","O","A","C"),"_se"))

df_bi <- merge(df_bi, m_theta.bfi_critA, by = "ID", all = TRUE)

(emp_rel_bfi_critA <- MFCblockInfo::calc.rel.emp(m_theta.bfi_critA[, 2:6], m_ses.bi_critA))
(emp_rel_bfi_critA_z <- psych::fisherz(sqrt(emp_rel_bfi_critA)))

(diff_bi_critA <- emp_rel_bfi_critA_z - emp_rel_bfi_all_z)

#--------------------------------- critB ####
m_theta.bfi_critB <- read.table("6_exploratory/fscores_bi_critB.dat", header = FALSE, na.strings = "*")[, c(71,61,63,65,67,69)]
colnames(m_theta.bfi_critB)<- c("ID", paste0(c("N","E","O","A","C"), "_bi_critB"))
m_ses.bi_critB  <- read.table("6_exploratory/fscores_bi_critB.dat", header = FALSE, na.strings = "*")[,c(62,64,66,68,70)]
colnames(m_ses.bi_critB )<- c(paste0(c("N","E","O","A","C"),"_se"))

df_bi <- merge(df_bi, m_theta.bfi_critB, by = "ID", all = TRUE)

(emp_rel_bfi_critB <- MFCblockInfo::calc.rel.emp(m_theta.bfi_critB[, 2:6], m_ses.bi_critB))
(emp_rel_bfi_critB_z <- psych::fisherz(sqrt(emp_rel_bfi_critB)))
  
(diff_bi_critB <-  emp_rel_bfi_critB_z - emp_rel_bfi_all_z)


#_________________________________ Relationship between Big Five, Narc with age & gender ______________________________ ####

explo_corr <- data.frame(var = colnames(df_bi[-c(1:3)]),
                         age_cor = NA, age_ci = NA, 
                         gender_cor = NA, gender_ci = NA, 
                         rel_est = NA, rel_est_z = NA)
explo_corr

# manifest ####
explo_corr$age_cor <- apply(df_bi[, 4:18], 2, FUN = function(x, a) numformat(round2(cor(x, a, use = "pairwise"),2)), a=df_bi$age)
explo_corr$age_ci <- apply(df_bi[, 4:18], 2, FUN = function(x, a) paste0("[", numformat(round2(cor.test(x, a)$conf.int[1], 2)),
                                                                        "; ", numformat(round2(cor.test(x, a)$conf.int[2], 2)), "]"), a = df_bi$age)
explo_corr$gender_cor <- apply(df_bi[, 4:18], 2, FUN = function(x, a) numformat(round2(cor(x, a, use = "pairwise"),2)), a = df_bi$gender)
explo_corr$gender_ci <- apply(df_bi[, 4:18], 2, FUN = function(x, a) paste0("[", numformat(round2(cor.test(x, a)$conf.int[1], 2)),
                                                                           "; ", numformat(round2(cor.test(x, a)$conf.int[2], 2)), "]"), a = df_bi$gender)
explo_corr$rel_est <- numformat(round2(c(emp_rel_bfi_all, emp_rel_bfi_critA, emp_rel_bfi_critB), 2))
explo_corr$rel_est_z <- numformat(round2(c(emp_rel_bfi_all_z, emp_rel_bfi_critA_z, emp_rel_bfi_critB_z),2))
explo_corr$rel_diff <- numformat(round2(c(rep(NA, 5), diff_bi_critA, diff_bi_critB), 2))

head(explo_corr)







# end 