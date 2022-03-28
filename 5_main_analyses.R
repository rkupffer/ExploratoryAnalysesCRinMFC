###################################################  Main Analyses ###################################################
df <- readRDS("4_indices.rds")

library(psych)
library(CRinMFC)
library(TirtAutomation)
library(tidyLPA)
library(apaTables)
library(MOTE)
library(mice)

# additional functions which are not part of any package
source("0_add_functions.R")

#####______________________________________  Descriptive Analyses ______________________________________ ####
indices <- c("time", "timeRTI", "SR_Effort", "SR_Attention", "SR_UseMe", "ir", "cs.mean", 
             "md.mean", "lcm.df", "lca.df", "sc.df", "tv.df", "missi")
indices.cutoff <- c("rt", "rti", "sr_e", "sr_a", "sr_u", "ir", "cs", "md", "lcm", "lca", "sc", "tv", "mis")
Index <- c("Response Time", "Response Time Index", "Self-Report Effort", "Self-Report Attention", "Self-Report UseMe", 
           "Instructed Response Triplet", "Consistency Score", "Mahalanobis D", "LongCombiMax", "LongCombiAvg",
           "SameCombi", "TripletVariance", "Missingness")
prop_table <- data.frame(Index, M = NA, SD = NA)
prop_table$M <- sapply(indices, FUN = function(indi) round2(mean(df[, indi], na.rm = TRUE), 2))
prop_table$SD <- sapply(indices, FUN = function(indi) round2(sd(df[, indi], na.rm = TRUE), 2))
# confidence interval
ci <- function(mean, sd, z){
  ci_up <- mean + z*sd
  ci_lo <- mean - z*sd
  paste0("[", round2(ci_lo, 2), "; ", round2(ci_up, 2), "]")
}
prop_table$CI <- ci(prop_table$M, prop_table$SD, z = 1.96)
prop_table$Cutoff <- c("mean-2sd", "mean+2sd", "3 - some", "3 - some", "1 - no",
                       "1", "mean-2sd", "2questD^2>Chi", "7", "2", ".4", ".7", "mean+2sd")
prop_table$percentage <- sapply(indices.cutoff, FUN = function(indi) round2(prop.table(table(df[, indi])), 2)[2])
prop_table[6, 6] <- .10 #percentage for at least one incorrectly responded irt
prop_table

#####______________________________________ RQ1: correlations and factor structure ______________________________________ ####
####_______ correlations ####
# with recoded indices (higher values = "more" careless responding) 
indices <- subset(df, select = c(time, timeRTI, SR_E_rec, SR_A_rec, SR_UseMe, ir, cs.mean_rec, md.mean, lcm.df, lca.df, sc.df, tv.df_rec, missi, tsc))
colnames(indices) <- c("time", "timeRTI", "SR_Eff", "SR_Att", "SR_Use", "irt", "cs", "md", "lcm", "lca", "sc", "tv", "missi", "tsc")
head(indices)

corr_table_indices <- cor(indices, use = "pairwise")
apa_corr <- apa(corr_table_indices , decimals = 2, leading = FALSE)
corr_table <- matrix(apa_corr, nrow = 14, ncol = 14)
corr_table[upper.tri(corr_table)] <- " "

####_______ factor structure ####
# preconditions
bartlett <- cortest.bartlett(corr_table_indices, n = nrow(df))
bartlett

kmo <- KMO(corr_table_indices)
round2(kmo$MSA, 2)
# preconditions for further analyses are not met


#####______________________________________  RQ2: Latent-Profile-Analysis ______________________________________ ####
# select indices we preregistered to use for lpa
indices_lpa <- subset(df, select = c(time, SR_Effort, SR_Attention, cs.mean, md.mean, lcm.df, lca.df, sc.df, tv.df))
colnames(indices_lpa) <- c("time","SR_E", "SR_A", "cs", "md", "lcm", "lca","sc", "tv")
head(indices_lpa)

# single imputation (m = 1) of missing values on variables: 
#  time (2%), Self-Report Effort(7%) and Attention(3%), and consistency score (3%)
indices_lpa_imp <- mice(indices_lpa, seed = 210705, maxit = 10, m = 1)
# plot(indices_lpa_imp)
indices_lpa_imp1 <- complete(indices_lpa_imp, action = 1)
# plausibility of values
range(indices_lpa_imp1[, "SR_E"])
range(indices_lpa_imp1[, "SR_A"])
range(indices_lpa_imp1[, "cs"])

# Latent profile analysis
lpa <- estimate_profiles(indices_lpa_imp1, 1:5, 
                      models = 2,
                      package = "mplus", 
                      ANALYSIS = "starts = 9000, 400;")
lpa

fit_table_lpa <- data.frame(get_fit(lpa)[,c("Model","Classes","BIC","Entropy","prob_min","prob_max")])
fit_table_lpa$Model <- 0:4
fit_table_lpa$BIC <- round2(fit_table_lpa$BIC, 0)
fit_table_lpa[,c("Entropy","prob_min","prob_max")] <- round2(fit_table_lpa[,c("Entropy","prob_min","prob_max")], 2)


lpa$model_2_class_2$warnings # NULL
lpa$model_2_class_3$warnings # NULL
lpa$model_2_class_4$warnings # NULL
lpa$model_2_class_5$warnings # "Less than 1% of cases were assigned to one of the profiles. 
# Interpret this solution with caution and consider other models."  
#"THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY. ESTIMATES CANNOT BE TRUSTED."

#estimates
est_lpa <- get_estimates(lpa)

est_mw <- subset(est_lpa, est_lpa$Category == "Means")
head(est_mw)

est_mod2.2 <- as.data.frame(lpa$model_2_class_2$dff)
est_mod2.3 <- as.data.frame(lpa$model_2_class_3$dff)
est_mod2.4 <- as.data.frame(lpa$model_2_class_4$dff)

v_index <- c(as.character(colnames(indices_lpa)), "missi", "rti", "SR_UseMe", "ir", "prop")

####_________   2-Class solution model_2_2   ####
# data frame with class means

mw_estim2.2 <- data.frame(matrix(NA, nrow = 2, ncol = length(v_index)))
colnames(mw_estim2.2) <- v_index
df$Class_mod2 <- est_mod2.2$Class

for (i in 1:nrow(mw_estim2.2)){
  mw_estim2.2[i, 1:9] <- apply(est_mod2.2[est_mod2.2$Class == i, v_index[1:9]], 2, mean, na.rm = TRUE)
  mw_estim2.2[i, 10:13] <- apply(df[df$Class_mod2 == i, c("missi", "rti", "SR_UseMe", "ir")], 2, mean, na.rm = TRUE)
  mw_estim2.2[i, 14] <- sum(est_mod2.2$Class == i)/nrow(df)
}
mw_estim2.2 <- round2(mw_estim2.2, 2)

#prepare for the paper: reorder the indices and export the table 
model1 <- t(mw_estim2.2[, c("prop","time","rti","SR_E","SR_A","SR_UseMe","cs","md","lcm","lca","sc","tv","ir","missi")])
colnames(model1) <- c("Class1", "Class2")
model1 

####_________   3-Class solution model_2_3   ####
mw_estim2.3 <- data.frame(matrix(NA, nrow = 3, ncol = length(v_index)))
colnames(mw_estim2.3) <- v_index
df$Class_mod3 <- est_mod2.3$Class

for (i in 1:nrow(mw_estim2.3)){
  mw_estim2.3[i, 1:9] <- apply(est_mod2.3[est_mod2.3$Class == i, v_index[1:9]], 2, mean, na.rm = TRUE)
  mw_estim2.3[i, 10:13] <- apply(df[df$Class_mod3 == i, c("missi", "rti", "SR_UseMe", "ir")], 2, mean, na.rm = TRUE)
  mw_estim2.3[i, 14] <- sum(est_mod2.3$Class == i)/nrow(df)
}
mw_estim2.3 <- round2(mw_estim2.3, 2)

#prepare for the paper: reorder the indices and export the table 
model2 <- t(mw_estim2.3[, c("prop","time","rti","SR_E","SR_A","SR_UseMe","cs","md","lcm","lca","sc","tv","ir","missi")])
model2 <- model2[, c(2,1,3)] #reorder: largest class should always be Class1
colnames(model2) <- c("Class1", "Class2", "Class3")
model2 


####_________   4-Class solution model_2_4  ####
mw_estim2.4 <- data.frame(matrix(NA, nrow = 4, ncol = length(v_index)))
colnames(mw_estim2.4) <- v_index
df$Class_mod4 <- est_mod2.4$Class

for (i in 1:nrow(mw_estim2.4)){
  mw_estim2.4[i, 1:9] <- apply(est_mod2.4[est_mod2.4$Class == i, v_index[1:9]], 2, mean, na.rm = TRUE)
  mw_estim2.4[i, 10:13] <- apply(df[df$Class_mod4 == i, c("missi", "rti", "SR_UseMe", "ir")], 2, mean, na.rm = TRUE)
  mw_estim2.4[i, 14] <- sum(est_mod2.4$Class == i)/nrow(df)
}
mw_estim2.4 <- round2(mw_estim2.4, 2)

#prepare for the paper: reorder the indices and classes and export the table 
model3 <- t(mw_estim2.4[, c("prop","time","rti","SR_E","SR_A","SR_UseMe","cs","md","lcm","lca","sc","tv","ir","missi")])
model3 <- model3[, c(2,1,4,3)]
colnames(model3) <- c("Class1", "Class2", "Class3", "Class4")
model3 

# save LPA results ####
save(indices_lpa, lpa,
     mw_estim2.2, mw_estim2.3, mw_estim2.4,
     file = "5_rq2_lpa_results.RData")


#####______________________________________ RQ3: proportion of CR for each index  ______________________________________ ####
# see section "Descriptive Analyses"

#####______________________________________ RQ4: Where does CR appear?  ______________________________________ ####
## 1) For both questionnaires (BFT and ORVIS) we will compute the proportion of careless respondents
#     identified by the indices
#     Consistency Score, Mahalanobis D, LongCombiMax, LongCombiAvg, SameCombi, and TripletVariance.
#     And as exploratory results: missing values and rti
indi <- c("cs", "md", "lcm", "lca", "sc", "tv", "rti", "mis")

ID <- df$ID
dat <- data.frame(ID)
head(dat,2)

### BFT ####
dat$cs.cf.bt <- ifelse(df$cs.bt < (mean(df$cs.bt, na.rm = TRUE) - 2*sd(df$cs.bt, na.rm = TRUE)), 1, 0)

dat$md.cf.bt <- df$md.bt.cr

dat$lcm.cf.bt <- ifelse(df$lcm.bt >= 7, 1, 0)
dat$lca.cf.bt <- ifelse(df$lca.bt >= 2, 1, 0)
dat$sc.cf.bt  <- ifelse(df$sc.bt >= .4, 1, 0)
dat$tv.cf.bt  <- ifelse(df$tv.bt <  .7, 1, 0)
dat$rti.cf.bt <- df$rti.bt
dat$mis.cf.bt <- ifelse(df$mis.bt > mean(df$mis.bt)+2*sd(df$mis.bt), 1, 0)

BFT <- round2(apply(dat[,-c(1)], 2, function(col) prop.table(table(col))[2]), 2)

### ORVIS ####
dat$cs.cf.os <- ifelse(df$cs.os < (mean(df$cs.os, na.rm = TRUE) - 2*sd(df$cs.os, na.rm = TRUE)), 1, 0)

dat$md.cf.os <- df$md.os.cr

dat$lcm.cf.os <- ifelse(df$lcm.os >= 7, 1, 0)
dat$lca.cf.os <- ifelse(df$lca.os >= 2, 1, 0)
dat$sc.cf.os  <- ifelse(df$sc.os  >= .4, 1, 0)
dat$tv.cf.os  <- ifelse(df$tv.os  <  .7, 1, 0)
dat$rti.cf.os <- df$rti.os
dat$mis.cf.os <- ifelse(df$mis.os > mean(df$mis.os)+2*sd(df$mis.os), 1, 0)

ORVIS <- round2(apply(dat[,grep("[.]os",colnames(dat))], 2, function(col) prop.table(table(col))[2]), 2)

Chi <- rep(NA, length(indi))
p <- rep(NA, length(indi))
rq4 <- data.frame(indi, BFT, ORVIS, Chi, p)
rq4

for(i in 1:length(indi)){
  conti_table <- table(dat[,paste0(indi[i],".cf.bt")], 
                 dat[,paste0(indi[i],".cf.os")])
  # expected frequencies must be equal or above 5
  # if not: continuity correction Yates
  concor <- ifelse(min(conti_table) < 5, TRUE, FALSE)
  mc <- mcnemar.test(conti_table, correct = concor)

  rq4[i,4] <- round2(mc$statistic,2)
  rq4[i,5] <- ifelse(mc$p.value<.001, "<.001", round2(mc$p.value,3))
}
rq4






#