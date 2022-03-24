########################                  RTI                                          #########################
#####                What is an adequate lower cut-off value for the response time index?                  #####
rm(list=ls())
library(CRinMFC)

# additional functions which are not part of any package
source("00_add_functions.R")


setwd("~/Documents/Git/CRinMFConline/z_otherScripts/RTI")
dat <- read.csv2("data_RTIinMFC_2021-03-19_18-58.csv", header = T, sep = ";", dec = ".", na.strings = "-9" )

dat <- dat[c(15:25),] #without old runs (tablets and phones included)
#case18 2021-02-26 15:25:57
#rerun:
#case19 2021-03-01 14:37:12
#same person, same computer, first run without double click, second time with double click
#remove rerun
dat <- dat[dat$STARTED != "2021-03-01 14:37:12", ]



# Section A - instruction: do not read the items, predefined rank order (1 2 3) ####
dat$mw1 <- rowMeans(dplyr::select(dat ,c(TIME004 : TIME014)))
MW_A <- mean(dat$mw1)
SD_A <- sd(dat$mw1)



# Section B - instruction: do not read the items, rank order not predefined ####
dat$mw2 <- rowMeans(dplyr::select(dat , c(TIME017 : TIME027)))
MW_B <- mean(dat$mw2) 
SD_B <- sd(dat$mw2)



# Section C - instruction: read the items, rank order not predefined ####
dat$mw3 <- rowMeans(dplyr::select(dat , c(TIME030: TIME040)))
MW_C <- mean(dat$mw3) 
SD_C <- sd(dat$mw3)


# Section D - instruction: read the items, rank the items according to how well they describe you####
dat$mw4 <- rowMeans(dplyr::select(dat , c(TIME043 :  TIME053)))
MW_D <- mean(dat$mw4) 
SD_D <- sd(dat$mw4)


# Data Frame erstellen ####
section <- c("A", "B", "C", "D")
read <- c("no", "no", "yes", "yes")
rank <- c("123", "equal", "equal","honestly")
df <- data.frame(section, read, rank, 
                 MW = round2(c(MW_A, MW_B, MW_C, MW_D),2), 
                 SD = round2(c(SD_A, SD_B, SD_C, SD_D),2))
df

t.test(dat$mw2, dat$mw3, paired = TRUE)   

# Figure ####
library(ggplot2)


#jpeg(file = "../z_otherScripts/RTI/ResponseTimeCondition.jpeg", width=430, height=330)

p <- ggplot(df, aes(x=section, y=MW)) + 
  geom_pointrange(aes(ymin=MW-SD, ymax=MW+SD))
print(p)
p + labs(y ="Average Time per Triplet (seconds)", x = "Condition")+ 
  theme_bw()

#dev.off()

              