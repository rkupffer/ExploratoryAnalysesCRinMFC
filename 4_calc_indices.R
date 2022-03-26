####________________________  calculation of the indices to detect careless responding in mfc questionnaires ________________________####
#df <- read.csv("DataExploratoryAnalysesCRinMFC.csv")[-1]
df <- readRDS("2_df_recoded.rds")

# which item belongs to which trait and how is it keyed
all_items <- read.csv2("2_All_Items_Coding.csv", sep = ",")[,2:8]

# packages
library(psych)
library(CRinMFC)
library(TirtAutomation)
library(MplusAutomation)
library(mice)

# additional functions which are not part of any package
source("0_add_functions.R")

# ------------------------------------------------- Response Time  -------------------------------------------------
time_vars <- df[,colnames(df[grepl("TIME", colnames(df))])]

# only time variables of pages during the questionnaires are needed
names_time_vars <- colnames(time_vars[,c(4:94,96:125)])

time.df <- df[,c(names_time_vars)]

avg_time_per_triplet <- rowMeans(time.df, na.rm = TRUE)
quants <- quantile(avg_time_per_triplet)

# replace values > 2*quartile by random uniform in interquartile range
count <- length(time.df[(time.df > 2*quants["75%"]) & (is.na(time.df)==FALSE)] )
set.seed(210629)
time.df[(time.df > 2*quants["75%"]) & (is.na(time.df)==FALSE)] <- runif(count, quants["25%"], quants["75%"])


df$time <- rowSums(time.df)/60

(100/(length(names_time_vars)*nrow(df)))*count # 4 percent of the values were replaced

# cutting value
df$rt <- ifelse(df$time < (mean(df$time, na.rm=TRUE) - 2*sd(df$time, na.rm=TRUE)), 1, 0)


####_______ Response Time Index --> RTI ####
timeRTI.df <- df[, c(names_time_vars)]

df$timeRTI <- responseTimeIndex(timeRTI.df, t_upper = 120, t_lower = 6)

# cutoff value
df$rti <- ifelse(df$timeRTI > (mean(df$timeRTI) + 2*sd(df$timeRTI)), 1, 0)
round2(prop.table(table(df$rti)),2)

# first and last questionnaire 
# bft
names_time_vars_bt <- df[, colnames(time_vars[,c(4:23)])]
df$timeRTI.bt <- responseTimeIndex(names_time_vars_bt, t_upper = 120, t_lower = 6)
df$rti.bt <- ifelse(df$timeRTI.bt > (mean(df$timeRTI.bt) + 2*sd(df$timeRTI.bt)), 1, 0)
round2(prop.table(table(df$rti.bt)),2)

# orvis
names_time_vars_os <- df[, colnames(time_vars[,c(96:125)])]
df$timeRTI.os <- responseTimeIndex(names_time_vars_os, t_upper = 120, t_lower = 6)
df$rti.os <- ifelse(df$timeRTI.os > (mean(df$timeRTI.os, na.rm = TRUE) + 2*sd(df$timeRTI.os, na.rm = TRUE)), 1, 0)
round2(prop.table(table(df$rti.os)),2)


# ------------------------------------------------- Self-Report Measures -------------------------------------------------
####_______ SR_Effort ####
# SR01_01 SRSI effort: I put forth ____ effort towards this study.
#1 - almost no  #2 - very little  #3 - some  #4 - quite a bit  #5 - a lot of
names(df)[names(df) == "SR01_01"] <- "SR_Effort"
table(df$SR_Effort)
df$sr_e <- ifelse(df$SR_Effort <= 3, 1, 0)
round2(prop.table(table(df$sr_e)), 2)

# recode (5 - almost no; 4 - very little; 3 - some; 2 - quite a bit; 1 - a lot of)
df$SR_E_rec <- 6 - df$SR_Effort
corr.test(df$SR_E_rec, df$SR_Effort)

####_______ SR_Attention ####
#SR02_01 SRSI attention: I gave this study ____ attention.
#1 - almost no  #2 - very little of my  #3 - some of my  #4 - most of my  #5 - my full
names(df)[names(df) == "SR02_01"] <- "SR_Attention"
table(df$SR_Attention)
df$sr_a <- ifelse(df$SR_Attention <= 3, 1, 0)
round2(prop.table(table(df$sr_a)), 2)

# recode
df$SR_A_rec <- 6 - df$SR_Attention
corr.test(df$SR_A_rec, df$SR_Attention)

####_______ SR_UseMe ####
#SR03: SRSI UseMe: Should your data be included in the analyses of the present study?
# 1 - yes  #2 - no
names(df)[names(df) == "SR03"] <- "SR_UseMe"
table(df$SR_UseMe)
df$SR_UseMe <- car::recode(df$SR_UseMe, recodes = "2 = 1; 1 = 0")
df$sr_u <- ifelse(df$SR_UseMe == 1, 1, 0)
round(prop.table(table(df$sr_u)), 2)


# ------------------------------------------------- Instructed Response Triplet  -------------------------------------------------
# 0 - correct answer; 1 - false answer; NA - no answer

#  VC01.1: Place this statement at the second rank.
#  VC01.2: Please place this statement at rank three.
#  VC01.3: Drag this statement to the top rank.
df$irt1 <- ifelse((is.na(df$VC01_01) & is.na(df$VC01_02) & is.na(df$VC01_03)), NA,
                    ifelse(df$VC01_01 ==2 & df$VC01_02 ==3 & df$VC01_03 ==1, 0, 1))
round2(prop.table(table(df$irt1)),2)

#  VC02.1: Please drag this box to the last rank.
#  VC02.2: Place this statement at rank one.
#  VC02.3: This statement should be placed at the second rank.
df$irt2 <- ifelse((is.na(df$VC02_01) & is.na(df$VC02_02) & is.na(df$VC02_03)), NA,
              ifelse(df$VC02_01 ==3 & df$VC02_02 ==1 & df$VC02_03 ==2, 0, 1))
round2(prop.table(table(df$irt2)),2)

df$ir <- rowSums(df[,which(names(df)=="irt1"):which(names(df)=="irt2")])

round2(prop.table(table(df$ir)),2)

# ------------------------------------------------- Consistency Score -------------------------------------------------

####_______ bft ---> bt ####
#b# blocks = 20
#m# items per block = 3

# import parameters from mplus output
bt <- readModels("3_Mplus/03_TIRT_BT.out", what=c("parameters", "sampstat"))
bt.pars <- bt$parameters$unstandardized

# factor scores
m_theta.bt <- read.table("3_Mplus/fscores_bt.dat", header = FALSE, na.strings = "*")[,c(71,61,63,65,67,69)]
colnames(m_theta.bt)<- c("ID",paste0(c("N","E","O","A","C"),"_bt"))
head(m_theta.bt, 3)

# keep IDs which are also included in df
m_theta.bt.all <- merge(m_theta.bt, df, by = "ID", all.y = TRUE)
m_theta.bt.ID <- subset(m_theta.bt.all, select = ID:C_bt)

m_theta.bt <- m_theta.bt.ID[, 2:6]

# prepare binary outcomes: only keep IDs which are also included in m_theta
d_bi.bt.sub <- subset(m_theta.bt.all, select = b.bt.1:b.bt.60)
dim(d_bi.bt.sub)

# Consistency Score: how many of the probabilities are above 0.5?
cs.bt <- consisScore(
  #mplus parameters
  quest.pars =   bt.pars,
  #number of blocks
  no.b = 20,
  #no.traits     - number of traits measured in the questionnaire
  no.traits = 5,
  #all.items     - item coding table
  all.items = all_items,
  #quest         - which questionnaire "BT"/"SD"/"BI"/"OS"...
  quest = "BT",
  #m.theta       - thetas
  m.theta = m_theta.bt,
  #d.bi          - coded df with binary outcomes (same length as m.fscores, hence sometimes a subsample)
  d.bi = d_bi.bt.sub)

m_theta.bt.ID$cs.bt <- cs.bt

# cs in df
df <- merge(df, m_theta.bt.ID[, c("ID","cs.bt")], by="ID", all.x = TRUE)

describe(df$cs.bt)
hist(df$cs.bt)

# recode
df$cs.bt_rec <- 1 - df$cs.bt



####_______ bfi ---> bi ####
#b# blocks = 20
#m# items per block = 3

# import parameters from mplus output
bi <- readModels("3_Mplus/03_TIRT_BI.out", what=c("parameters","sampstat"))
bi.pars <- bi$parameters$unstandardized

# factor scores
m_theta.bi <- read.table("3_Mplus/fscores_bi.dat", header = FALSE, na.strings = "*")[,c(71,61,63,65,67,69)]
colnames(m_theta.bi)<- c("ID",paste0(c("N","E","O","A","C"),"_bi"))

# keep IDs which are also included in df
m_theta.bi.all <- merge(m_theta.bi, df, by = "ID", all.y = TRUE)
m_theta.bi.ID <- subset(m_theta.bi.all, select = ID:C_bi)
head(m_theta.bi.ID,3)
dim(m_theta.bi.ID)

m_theta.bi <- m_theta.bi.ID[,2:6]

# prepare binary outcomes: only keep IDs which are also included in m_theta
d_bi.bi.sub <- subset(m_theta.bi.all, select = b.bi.1:b.bi.60)

# Consistency Score: how many of the probabilities are above 0.5?
cs.bi <- consisScore(
  #mplus parameters
  quest.pars =   bi.pars,
  #number of blocks
  no.b = 20,
  #no.traits     - number of traits measured in the questionnaire
  no.traits = 5,
  #all.items     - item coding table
  all.items = all_items,
  #quest         - which questionnaire "BT"/"SD"/"BI"/"OS"...
  quest = "BI",
  #m.theta       - thetas
  m.theta = m_theta.bi,
  #d.bi          - coded df with binary outcomes (same length as m.fscores, hence sometimes a subsample)
  d.bi = d_bi.bi.sub)

m_theta.bi.ID$cs.bi <- cs.bi

#cs in df
df <- merge(df, m_theta.bi.ID[,c("ID","cs.bi")], by="ID", all.x = TRUE)

describe(df$cs.bi)
hist(df$cs.bi)

df$cs.bi_rec <- 1 - df$cs.bi


####_______ HEXACO ---> ho ####
#b# blocks = 20

# import parameters from mplus output
ho <- readModels("3_Mplus/03_TIRT_HO.out", what=c("parameters","sampstat"))
ho.pars <- ho$parameters$unstandardized

# factor scores
m_theta.ho <- read.table("3_Mplus/fscores_ho.dat", header = FALSE, na.strings = "*")[,c(73,61,63,65,67,69,71)]
colnames(m_theta.ho)<- c("ID",paste0(c("H","E","X","A","C","O"),"_ho"))
dim(m_theta.ho)

# keep IDs which are also included in df
m_theta.ho.all <- merge(m_theta.ho, df, by = "ID", all.y = TRUE)
m_theta.ho.ID <- subset(m_theta.ho.all, select = ID:O_ho)
head(m_theta.ho.ID,3)
dim(m_theta.ho.ID)

m_theta.ho <- m_theta.ho.ID[,2:7]

# prepare binary outcomes: only keep IDs which are also included in m_theta
d_bi.ho.sub <- subset(m_theta.ho.all, select = b.ho.1:b.ho.60)

# Consistency Score: how many of the probabilities are above 0.5?
cs.ho <- consisScore(
  #mplus parameters
  quest.pars =   ho.pars,
  #number of blocks
  no.b = 20,
  #no.traits     - number of traits measured in the questionnaire
  no.traits = 6,
  #all.items     - item coding table
  all.items = all_items,
  #quest         - which questionnaire "BT"/"SD"/"BI"/"OS"...
  quest = "HO",
  #m.theta       - thetas
  m.theta = m_theta.ho,
  #d.bi          - coded df with binary outcomes (same length as m.fscores, hence sometimes a subsample)
  d.bi = d_bi.ho.sub)

m_theta.ho.ID$cs.ho <- cs.ho

# cs in df
df <- merge(df, m_theta.ho.ID[,c("ID","cs.ho")], by="ID", all.x = TRUE)

describe(df$cs.ho)
hist(df$cs.ho)

df$cs.ho_rec <- 1 - df$cs.ho

####_______ IPIP ---> ip ####
#b# blocks = 20

# import parameters from mplus output
ip <- readModels("3_Mplus/03_TIRT_IP.out", what=c("parameters","sampstat"))
ip.pars <- ip$parameters$unstandardized

# factor scores
m_theta.ip <- read.table("3_Mplus/fscores_ip.dat", header = FALSE, na.strings = "*")[,c(73,61,63,65,67,69,71)]
colnames(m_theta.ip)<- c("ID",paste0(c("C","D","I","J","R","S"),"_ip"))
# keep IDs which are also included in df
m_theta.ip.all <- merge(m_theta.ip, df, by = "ID", all.y = TRUE)
m_theta.ip.ID <- subset(m_theta.ip.all, select = ID:S_ip)
head(m_theta.ip.ID,3)
dim(m_theta.ip.ID)

m_theta.ip <- m_theta.ip.ID[,2:7]

# prepare binary outcomes: only keep IDs which are also included in m_theta
d_bi.ip.sub <- subset(m_theta.ip.all, select = b.ip.1:b.ip.60)


# Consistency Score: how many of the probabilities are above 0.5?
cs.ip <- consisScore(
  #mplus parameters
  quest.pars =   ip.pars,
  #number of blocks
  no.b = 20,
  #no.traits     - number of traits measured in the questionnaire
  no.traits = 6,
  #all.items     - item coding table
  all.items = all_items,
  #quest         - which questionnaire "BT"/"SD"/"BI"/"OS"...
  quest = "IP",
  #m.theta       - thetas
  m.theta = m_theta.ip,
  #d.bi          - coded df with binary outcomes (same length as m.fscores, hence sometimes a subsample)
  d.bi = d_bi.ip.sub)

m_theta.ip.ID$cs.ip <- cs.ip

# cs in df
df <- merge(df, m_theta.ip.ID[,c("ID","cs.ip")], by="ID", all.x = TRUE)

describe(df$cs.ip)
hist(df$cs.ip)

df$cs.ip_rec <- 1 - df$cs.ip


####_______ ORVIS ---> os ####
#b# blocks = 29 
#one triplet was deleted due to convergence problems of the thurstonian irt model

# import parameters from mplus output
os <- readModels("3_Mplus/03_TIRT_OS.out", what=c("parameters","sampstat"))
os.pars <- os$parameters$unstandardized

#factor scores
m_theta.os <- read.table("3_Mplus/fscores_os.dat", header = FALSE, na.strings = "*")[,c(104,88,90,92,94,96,98,100,102)]
colnames(m_theta.os)<- c("ID",paste0(c("L","O","A","C","N","P","D","E"),"_os"))
head(m_theta.os,3)
# keep IDs which are also included in df
m_theta.os.all <- merge(m_theta.os, df, by = "ID", all.y = TRUE)
m_theta.os.ID <- subset(m_theta.os.all, select = ID:E_os)
head(m_theta.os.ID,3)
dim(m_theta.os.ID)

m_theta.os <- m_theta.os.ID[,2:9]

# prepare binary outcomes: only keep IDs which are also included in m_theta
# without triplet 12
d_bi.os.sub2 <- subset(m_theta.os.all, select = c(b.os.1:b.os.87)) 

# remove deleted triplet form all_items
all_items.os <- all_items[which(all_items$questionnaire == "OS" & all_items$triplet != 12),]

# Consistency Score: how many of the probabilities are above 0.5?
cs.os <- consisScore(
  #mplus parameters
  quest.pars =   os.pars,
  #number of blocks
  no.b = 29,
  #no.traits     - number of traits measured in the questionnaire
  no.traits = 8,
  #all.items     - item coding table
  all.items = all_items.os,
  #quest         - which questionnaire "BT"/"SD"/"BI"/"OS"...
  quest = "OS",
  #m.theta       - thetas
  m.theta = m_theta.os,
  #d.bi          - coded df with binary outcomes (same length as m.fscores, hence sometimes a subsample)
  d.bi = d_bi.os.sub2)

m_theta.os.ID$cs.os <- cs.os

#cs in df
df <- merge(df, m_theta.os.ID[,c("ID","cs.os")], by="ID", all.x = TRUE)

describe(df$cs.os)
hist(df$cs.os)

df$cs.os_rec <- 1 - df$cs.os


####_______ rowMeansCS ####
df$cs.mean <- rowMeans(df[, c("cs.bt","cs.bi","cs.ho","cs.ip","cs.os")], na.rm = TRUE)

df$cs.mean_rec <- rowMeans(df[, c("cs.bt_rec","cs.bi_rec","cs.ho_rec","cs.ip_rec","cs.os_rec")], na.rm = TRUE)
describe(df$cs.mean)
hist(df$cs.mean)

df$cs <- ifelse(df$cs.mean < (mean(df$cs.mean, na.rm=TRUE) - 2*sd(df$cs.mean, na.rm=TRUE)), 1, 0)

round2(prop.table(table(df$cs)),2)

# ------------------------------------------------- MahalanobisD -------------------------------------------------
load("4_imputed.Rdata") # to not have to run the imputations every single time
#lines are commented out in order to avoid running them again
####_______ bft ---> bt ####
cor.bt <- cor2mat(bt)
thresh.bt <- thresh2vec(bt)
d_center.bt <- sweep(d_bi.bt.sub, 2L, thresh.bt)
# impute missing values using mice 
# d_center.bt.imp <- mice(d_center.bt, seed = 210702, m = 10, maxit = 10, method = "pmm")
# did the imputation converged? 
# "little trend and the streams mingle very well right from the start" (Buuren & Groothuis-Oudshoorn, 2010, p.39)
# plot(d_center.bt.imp) #(some plots cause R 4.1.1 on Linux to crash)
# calculate md in each of the imputed data sets
md.solutions.bt <- matrix(NA, nrow = nrow(d_center.bt), ncol = 10)
for(i in 1:10){
  complete.data.no.X <- mice::complete(d_center.bt.imp, action = i)
  md.solutions.bt[,i] <- mahalanobis(complete.data.no.X, center = FALSE, cov = cor.bt)
}
head(md.solutions.bt)
df$md.bt <- rowMeans(md.solutions.bt)
describe(df$md.bt)
df$md.bt.cr <- ifelse(abs(df$md.bt) > qchisq(p = 1-0.05, df = 60), 1, 0)


####_______ bfi ---> bi ####
cor.bi <- cor2mat(bi)
thresh.bi <- thresh2vec(bi)
d_center.bi <- sweep(d_bi.bi.sub, 2L, thresh.bi)
# impute missing values using mice 
# d_center.bi.imp <- mice(d_center.bi, seed = 210702, m = 10, maxit = 10, method = "pmm")
# did the imputation converged? 
# "little trend and the streams mingle very well right from the start" (Buuren & Groothuis-Oudshoorn, 2010, p.39)
# plot(d_center.bi.imp)
# calculate md in each of the imputed data sets
md.solutions.bi <- matrix(NA, nrow = nrow(d_center.bi), ncol = 10)
for(i in 1:10){
  complete.data.no.X <- mice::complete(d_center.bi.imp, action = i)
  md.solutions.bi[,i] <- mahalanobis(complete.data.no.X, center = FALSE, cov = cor.bi)
}
head(md.solutions.bi)
df$md.bi <- rowMeans(md.solutions.bi)
describe(df$md.bi)
df$md.bi.cr <- ifelse(abs(df$md.bi) > qchisq(p = 1-0.05, df = 60), 1, 0)

####_______ HEXACO ---> ho ####
cor.ho <- cor2mat(ho)
thresh.ho <- thresh2vec(ho)
d_center.ho <- sweep(d_bi.ho.sub, 2L, thresh.ho)
# impute missing values using mice 
# d_center.ho.imp <- mice(d_center.ho, seed = 210702, m = 10, maxit = 10, method = "pmm")
# did the imputation converged? 
# "little trend and the streams mingle very well right from the start" (Buuren & Groothuis-Oudshoorn, 2010, p.39)
# plot(d_center.ho.imp)
# calculate md in each of the imputed data sets
md.solutions.ho <- matrix(NA, nrow = nrow(d_center.ho), ncol = 10)
for(i in 1:10){
  complete.data.no.X <- mice::complete(d_center.ho.imp, action = i)
  md.solutions.ho[,i] <- mahalanobis(complete.data.no.X, center = FALSE, cov = cor.ho)
}
head(md.solutions.ho)
df$md.ho <- rowMeans(md.solutions.ho)
describe(df$md.ho)
df$md.ho.cr <- ifelse(abs(df$md.ho) > qchisq(p = 1-0.05, df = 60), 1, 0)

####_______ IPIP ---> ip ####
cor.ip <- cor2mat(ip)
thresh.ip <- thresh2vec(ip)
d_center.ip <- sweep(d_bi.ip.sub, 2L, thresh.ip)
# impute missing values using mice 
# d_center.ip.imp <- mice(d_center.ip, seed = 210702, m = 10, maxit = 10, method = "pmm")
# did the imputation converged? 
# "little trend and the streams mingle very well right from the start" (Buuren & Groothuis-Oudshoorn, 2010, p.39)
# plot(d_center.ip.imp)
# calculate md in each of the imputed data sets
md.solutions.ip <- matrix(NA, nrow = nrow(d_center.ip), ncol = 10)
for(i in 1:10){
  complete.data.no.X <- mice::complete(d_center.ip.imp, action = i)
  md.solutions.ip[,i] <- mahalanobis(complete.data.no.X, center = FALSE, cov = cor.ip)
}
head(md.solutions.ip)
df$md.ip <- rowMeans(md.solutions.ip)
describe(df$md.ip)
df$md.ip.cr <- ifelse(abs(df$md.ip) > qchisq(p = 1-0.05, df = 60), 1, 0)

####_______ ORVIS ---> os ####
cor.os <- cor2mat(os)
thresh.os <- thresh2vec(os)
d_center.os <- sweep(d_bi.os.sub2, 2L, thresh.os)
# impute missing values using mice 
# d_center.os.imp <- mice(d_center.os, seed = 210702, m = 10, maxit = 10, method = "pmm")
# did the imputation converged? 
# "little trend and the streams mingle very well right from the start" (Buuren & Groothuis-Oudshoorn, 2010, p.39)
# plot(d_center.os.imp)
# calculate md in each of the imputed data sets
md.solutions.os <- matrix(NA, nrow = nrow(d_center.os), ncol = 10)
for(i in 1:10){
  complete.data.no.X <- mice::complete(d_center.os.imp, action = i)
  md.solutions.os[,i] <- mahalanobis(complete.data.no.X, center = FALSE, cov = cor.os)
}
head(md.solutions.os)
df$md.os <- rowMeans(md.solutions.os)
describe(df$md.os)
df$md.os.cr <- ifelse(abs(df$md.os) > qchisq(p = 1-0.05, df = 87), 1, 0)

####_______ rowSumsMD ####
df$md.sum <- rowSums(df[, c("md.bt.cr","md.bi.cr","md.ho.cr","md.ip.cr","md.os.cr")])
df$md.mean <- rowMeans(df[, c("md.bt","md.bi","md.ho","md.ip","md.os")])
table(df$md.sum )
df$md <- ifelse(df$md.sum >= 2, 1, 0)
round2(prop.table(table(df$md)),2)

# ------------------------------------------------- Long- and sameOrder -------------------------------------------------
####_______ bft ---> bt ####
# create blocks of possible rank orders
d_b_o.bt <- sortOrders(no.b = 20, no.m = 3,
                d.r = subset(df, select = MT01_01:MT20_03),  #d.r - data.frame with rawdata (unrecoded answers to triplets)
                d.k = subset(df, select = KE01_01:KE01_20))  #d.k - data.frame with key variables
head(d_b_o.bt)

# calculate all longOrders (and save the results in a matrix)
m_lcm.bt <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.bt,     #df with ordered blocks of rank orders
                      no.b = 20)          #number of blocks
head(m_lcm.bt)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.bt <- longOrderMax(m_lcm.bt)
table(df$lcm.bt)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.bt <- longOrderAvg(m_lcm.bt)
hist(df$lca.bt,20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.bt <- sameOrder(d_b_o.bt, 20)
hist(df$sc.bt)

####_______ bfi ---> bi ####
# create blocks of possible rank orders
d_b_o.bi <- sortOrders(no.b = 20, no.m = 3,
                   d.r = subset(df, select = BF01_01:BF20_03),  #d.r - data.frame with rawdata (unrecoded answers to triplets)
                   d.k = subset(df, select = KE06_01:KE06_20))  #d.k - data.frame with key variables
head(d_b_o.bi)

# calculate all longOrders (and save the results in a matrix)
m_lcm.bi <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.bi,     #df with ordered blocks of rank orders
                      no.b = 20)            #number of blocks
head(m_lcm.bi)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.bi <- longOrderMax(m_lcm.bi)
table(df$lcm.bi)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.bi <- longOrderAvg(m_lcm.bi)
hist(df$lca.bi, 20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.bi <- sameOrder(d_b_o.bi, 20)
hist(df$sc.bi)

####_______ sd3 ---> sd ####
# create blocks of possible rank orders
d_b_o.sd <- sortOrders(no.b = 10, # an IRT appeared at a random position within the blocks --> 10 blocks
                   no.m = 3,
                   d.r = subset(df, select = c(DT01_01:DT09_03,
                                               VC01_01:VC01_03)),  #IRT1
                   d.k = subset(df, select = KE02_01:KE02_10))

head(d_b_o.sd)

# calculate all longOrders (and save the results in a matrix)
m_lcm.sd <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.sd,     #df with ordered blocks of rank orders
                      no.b = 10)            #number of blocks
head(m_lcm.sd)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.sd <- longOrderMax(m_lcm.sd)
table(df$lcm.sd)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.sd <- longOrderAvg(m_lcm.sd) 
hist(df$lca.sd, 20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.sd <- sameOrder(d_b_o.sd, 10)
hist(df$sc.sd)

####_______ HEXACO ---> ho ####
# create blocks of possible rank orders
d_b_o.ho <- sortOrders(no.b = 20, no.m = 3,
                   d.r = subset(df, select = HM01_01:HM20_03),  #d.r - data.frame with rawdata (unrecoded answers to triplets)
                   d.k = subset(df, select = KE03_01:KE03_20))  #d.k - data.frame with key variables
head(d_b_o.ho)

# calculate all longOrders (and save the results in a matrix)
m_lcm.ho <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.ho,     #df with ordered blocks of rank orders
                      no.b = 20)            #number of blocks
head(m_lcm.ho)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.ho <- longOrderMax(m_lcm.ho)
table(df$lcm.ho)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.ho <- longOrderAvg(m_lcm.ho)
hist(df$lca.ho, 20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.ho <- sameOrder(d_b_o.ho, 20)
hist(df$sc.ho)

####_______ IPIP ---> ip ####

# create blocks of possible rank orders
d_b_o.ip <- sortOrders(no.b = 21, # an IRT appeared at a random position within the blocks --> 21 blocks in total
                   no.m = 3,
                   d.r = subset(df, select = c(IP01_01:IP20_03, 
                                               VC02_01:VC02_03)),  # IRT2
                   d.k = subset(df, select = KE04_01:KE04_21))
head(d_b_o.ip)

# calculate all longOrders (and save the results in a matrix)
m_lcm.ip <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.ip,     #df with ordered blocks of rank orders
                      no.b = 21)          #number of blocks
head(m_lcm.ip)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.ip <- longOrderMax(m_lcm.ip)
table(df$lcm.ip)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.ip <- longOrderAvg(m_lcm.ip)
hist(df$lca.ip, 20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.ip <- sameOrder(d_b_o.ip, 21)
hist(df$sc.ip)


####_______ ORVIS ---> os ####
# create blocks of possible rank orders
d_b_o.os <- sortOrders(no.b = 30, no.m = 3,
                   d.r = subset(df, select = OR01_01:OR30_03),  #d.r - data.frame with rawdata (unrecoded answers to triplets)
                   d.k = subset(df, select = KE05_01:KE05_30))  #d.k - data.frame with key variables
head(d_b_o.os)

# calculate all longOrders (and save the results in a matrix)
m_lcm.os <- longOrder(IDs = df[, "ID"],     #df with ID variable
                      d.b.o = d_b_o.os,     #df with ordered blocks of rank orders
                      no.b = 30)          #number of blocks
head(m_lcm.os)

# longOrderMax: calculate the maximum of the six longOrder values for every person
df$lcm.os <- longOrderMax(m_lcm.os)
table(df$lcm.os)

# longOrderAvg: calculate the average over all six longOrder values for every person
df$lca.os <- longOrderAvg(m_lcm.os)
hist(df$lca.os, 20)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.os <- sameOrder(d_b_o.os, 30)
hist(df$sc.os)

####_______ whole survey ####
#### ordered df
# create blocks of possible rank orders
df.raw <- subset(df, select = c(MT01_01:MT20_03, #d.r.bt
                                BF01_01:BF20_03, #d.r.bi
                                DT01_01:DT09_03, #d.r.sd
                                VC01_01:VC01_03, #d.r.irt1
                                HM01_01:HM20_03, #d.r.ho
                                IP01_01:IP20_03, #d.r.ip
                                VC02_01:VC02_03, #d.r.irt2
                                OR01_01:OR30_03))#d.r.os

dim(df.raw)# 121 blocks * 3 items = 363 check

df.key <- subset(df, select = c(KE01_01:KE01_20, #d.k.bt
                                KE06_01:KE06_20, #d.k.bi
                                KE02_01:KE02_10, #d.k.sd
                                KE03_01:KE03_20, #d.k.ho
                                KE04_01:KE04_21, #d.k.ip
                                KE05_01:KE05_30))#d.k.os

dim(df.key)# 121 blocks --> check

d_b_o.df <- sortOrders(no.b = 121, no.m = 3, d.r = df.raw, d.k = df.key)
head(d_b_o.df, 3)

# calculate all longOrders (and save the results in a matrix)
m_lcm.df <- longOrder(IDs = df[, "ID"], 
                      d.b.o = d_b_o.df, 
                      no.b = 121)
head(m_lcm.df)


# longOrderMax: maximum of the six longOrder values for every person
df$lcm.df <- longOrderMax(m_lcm.df)
table(df$lcm.df)
# cutoff
df$lcm <- ifelse(df$lcm.df >= 7, 1, 0)

round2(prop.table(table(df$lcm)),2)


# longOrderAvg: average over all six longOrder values for every person
df$lca.df <- longOrderAvg(m_lcm.df)
hist(df$lca.df)
table(df$lca.df)

# cutoff
df$lca <- ifelse(df$lca.df >= 2, 1, 0)

round2(prop.table(table(df$lca)),2)

# sameOrder: proportion of triplets for which a person copied the presented order
df$sc.df <- sameOrder(d_b_o.df, 121)
hist(df$sc.df)

#cutoff
df$sc <- ifelse(df$sc.df >= .4, 1, 0)

round2(prop.table(table(df$sc)),2)



# which order chosen the most often?
round2(mean(m_lcm.df$c123),2)#2.06 --> 1 2 3 (on average 2 times in a row)
round2(mean(m_lcm.df$c132),2)#1.39
round2(mean(m_lcm.df$c213),2)#1.58
round2(mean(m_lcm.df$c231),2)#1.03
round2(mean(m_lcm.df$c312),2)#1.33
round2(mean(m_lcm.df$c321),2)#1.04

# ------------------------------------------------- tripletVariance -------------------------------------------------
####_______ bft ---> bt ####
df$tv.bt <- tripletVariance(d_b_o.bt)
hist(df$tv.bt)

# recode
df$tv.bt_rec <- 1 - df$tv.bt
hist(df$tv.bt_rec)

####_______ bfi ---> bi ####
df$tv.bi <- tripletVariance(d_b_o.bi)
hist(df$tv.bi)

df$tv.bi_rec <- 1 - df$tv.bi

####_______ sd3 ---> sd ####
df$tv.sd <- tripletVariance(d_b_o.sd)
hist(df$tv.sd)

df$tv.sd_rec <- 1 - df$tv.sd

####_______ HEXACO ---> ho ####
df$tv.ho <- tripletVariance(d_b_o.ho)
hist(df$tv.ho)

df$tv.ho_rec <- 1 - df$tv.ho

####_______ IPIP ---> ip ####
df$tv.ip <- tripletVariance(d_b_o.ip)
hist(df$tv.ip)

df$tv.ip_rec <- 1 - df$tv.ip

####_______ ORVIS ---> os ####
df$tv.os <- tripletVariance(d_b_o.os)
hist(df$tv.os)

df$tv.os_rec <- 1 - df$tv.os

####_______ whole survey ####
df$tv.df <- tripletVariance(d_b_o.df)
describe(df$tv.df)
hist(df$tv.df)

df$tv.df_rec <- 1 - df$tv.df

# cutoff
df$tv <- ifelse(df$tv.df < .7, 1, 0)

round2(prop.table(table(df$tv)), 2)


# ------------------------------------------------- TotalSumScore -------------------------------------------------
df$tsc <- rowSums(df[, c("rt","sr_e","sr_a","sr_u","ir","cs","md","lcm","lca","sc","tv")], na.rm = TRUE)
table(df$tsc)


# ------------------------------------------------- missing values  -------------------------------------------------
# counts how many missing values a person has on all questionnaire pages
df$missi <- rowMeans(is.na(df.raw))
hist(df$missi)
df$mis.bt <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "MT01_01"): which(colnames(df.raw) == "MT20_03")]))
df$mis.bi <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "BF01_01"): which(colnames(df.raw) == "BF20_03")]))
df$mis.sd <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "DT01_01"): which(colnames(df.raw) == "DT09_03")]))
df$mis.ho <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "HM01_01"): which(colnames(df.raw) == "HM20_03")]))
df$mis.ip <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "IP01_01"): which(colnames(df.raw) == "IP20_03")]))
df$mis.os <- rowMeans(is.na(df.raw[, which(colnames(df.raw) == "OR01_01"): which(colnames(df.raw) == "OR30_03")]))
df$mis.ir <- rowMeans(is.na(df.raw[, c(which(colnames(df.raw) == "VC01_01"): which(colnames(df.raw) == "VC01_03"),
                                       which(colnames(df.raw) == "VC02_01"): which(colnames(df.raw) == "VC02_03"))]))
df$mis <- ifelse(df$missi > mean(df$missi)+2*sd(df$missi), 1, 0)

# ------------------------------------------------- dataframe with ID and CR indices  -------------------------------------------------

indi <- df[,c("ID",
              "time", "rt", "timeRTI", "rti", "rti.bt", "rti.os",
              "SR_Effort", "SR_E_rec", "sr_e",
              "SR_Attention","SR_A_rec", "sr_a",
              "SR_UseMe", "sr_u",
              "irt1", "irt2", "ir", 
              "cs.bt", "cs.bi", "cs.ho", "cs.ip", "cs.os", "cs.mean", "cs.mean_rec", "cs",
              "md.bt", "md.bi", "md.ho", "md.ip", "md.os", 
              "md.bt.cr","md.bi.cr","md.ho.cr","md.ip.cr","md.os.cr", "md.sum", "md.mean", "md",
              "lcm.bt", "lcm.bi", "lcm.sd", "lcm.ho", "lcm.ip", "lcm.os", "lcm.df", "lcm",
              "lca.bt", "lca.bi", "lca.sd", "lca.ho", "lca.ip", "lca.os", "lca.df", "lca",
              "sc.bt", "sc.bi", "sc.sd", "sc.ho", "sc.ip", "sc.os", "sc.df", "sc",
              "tv.bt", "tv.bi", "tv.sd", "tv.ho", "tv.ip", "tv.os", "tv.df", "tv", "tv.df_rec",
              "mis.bt", "mis.bi", "mis.sd", "mis.ho", "mis.ip", "mis.os", "mis.ir", "mis", "missi",
              "tsc")]
dim(indi)
head(indi, 3)

#### save
saveRDS(indi, file = "4_indices.rds")

