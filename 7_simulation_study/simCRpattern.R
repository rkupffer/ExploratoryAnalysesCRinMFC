#### simulation and detection of different careless responding pattern  ####

# In this script different careless and careful responding pattern 
# will be simulated with varying proportions of careless respondents and 
# varying proportions of different response pattern.

# packages ----
library(MFCblockInfo)
library(TirtAutomation)
library(CRinMFC)
library(doParallel)

rm(list = ls())

# function: string ranks to numeric ----
ranks2num<- function(rankOrder){
  out <- rbind(as.numeric(substr(rankOrder, 1, 1)),
               as.numeric(substr(rankOrder, 2, 2)),
               as.numeric(substr(rankOrder, 3, 3)))
}

# function: sample and repeat rank orders ----

sampRepRank <- function(ranks2sample,    # list with rank orders as strings
                        m_out,           # output matrix
                        b,               # no of blocks
                        no_ranks2sample, # no of ranks drawn (per participant)
                        no_rep){         # number of repetitions per rank
  
  # strong repetition of rank orders - stroRepOrder  ------------------------
  if(no_ranks2sample == no_rep){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(ranks2sample, MARGIN = 2, FUN = sample,
                        size = nrow(m_out), replace = TRUE)
    
    # character to numeric ranks
    e1 <- as.numeric(substr(rankOrders, 1, 1))
    e2 <- as.numeric(substr(rankOrders, 2, 2))
    e3 <- as.numeric(substr(rankOrders, 3, 3))
    
    # repeat rank order no_rep times
    m_out <- rep(cbind(e1, e2, e3), no_rep)
  }
  
  # moderate repetition of rank orders - modRepOrder ------------------------
  if(no_ranks2sample == b/no_rep){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(rankOrders, 
                        MARGIN = 1:2, 
                        function(x) sample(ranks2sample, 
                                           size = 1, 
                                           replace = TRUE))
    
    
    # character to numeric ranks
    e1 <- matrix(data = NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    e1[ , 1:no_ranks2sample] <- as.numeric(
      substr(rankOrders[ , 1:no_ranks2sample], 1, 1))
    
    e2 <- matrix(NA, nrow(m_out), no_ranks2sample)
    e2[ , 1:no_ranks2sample] <- as.numeric(
      substr(rankOrders[ , 1:no_ranks2sample], 2, 2))
    
    e3 <- matrix(NA, nrow(m_out), no_ranks2sample)
    e3[ , 1:no_ranks2sample]  <- as.numeric(
      substr(rankOrders[ , 1:no_ranks2sample], 3, 3))
    
    # repeat each of the rank orders no_rep times
    m_ro1 <- cbind(e1[ , 1], e2[ , 1], e3[ , 1])
    m_ro2 <- cbind(e1[ , 2], e2[ , 2], e3[ , 2])
    m_ro3 <- cbind(e1[ , 3], e2[ , 3], e3[ , 3])
    m_ro4 <- cbind(e1[ , 4], e2[ , 4], e3[ , 4])
    
    m_out <- matrix(data = cbind(replicate(no_rep, m_ro1),
                                 replicate(no_rep, m_ro2),
                                 replicate(no_rep, m_ro3),
                                 replicate(no_rep, m_ro4)),
                    ncol = ncol(m_out), nrow = nrow(m_out),
                    byrow = FALSE)  
  }
  
  # random selection of rank orders - randomOrder  ------------------------
  if(no_ranks2sample == b){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(rankOrders, 
                        MARGIN = 1, 
                        function(x) sample(ranks2sample, 
                                           size = b, 
                                           replace = TRUE))

    # character to numeric ranks
    m_out_num <- matrix(data = NA, nrow = nrow(m_out), ncol = no_ranks2sample*3)
    
    m_out_num <- apply(rankOrders,
                       MARGIN = 1,
                       FUN = ranks2num)
    
    m_out <- t(m_out_num)

  }
}

####------------------- simulation design -------------------------####

# fixed parameter ----

# sample and BFT questionnaire meta data ----
all_items <- read.csv2("../ExploratoryAnalysesCRinMFC/2_All_Items_Coding.csv", 
                       sep = ",")[,2:8]

# number of blocks in the questionnaire
b <- 20

# number of items per block
m <- 3

# sample size
N <- 1000

# rank orders stored in a matrix
ro <- matrix(c("123", "132", "213", "231", "312", "321"), 6, byrow = TRUE)

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "BT", no.traits = 5)

# variance-covariance matrix (Anglim et al., 2020, p. 60)
m_sigma <- matrix(data = c(1, -.49, -.19, -.17, -.48, 
                           -.49, 1, .30, .08, .19, 
                           -.19, .30, 1, .20, .14,
                           -.17, .08, .20, 1, .32,
                           -.48, .19, .14, .32, 1), 
                  nrow = 5, ncol = 5, byrow = TRUE)

# vector with means fixed to 0
v_mu <- rep(0, 5)

# range of factor loadings
load <- c(-1, 1)

# range of intercepts
int <- c(.65, .95)

# varying parameter ----

# test design
p <- seq(0, 1, .25)
factor_randomOrder <- p
factor_strongRepOrder <- p
factor_moderateRepOrder <- p
factor_propOfCR <- seq(.02, .4, by = .05)

# smaller number of replications for test
r <- 1:10

# smaller sample size for test
N <- 10




# actual design
#factor_randomOrder <- c(rep(c(0,.25,.5,.75,1),2))
#factor_strongRepOrder <- c(1,.75,.5,.25, rep(0, 6))
#factor_moderateRepOrder <- c(rep(0, 5), 1,.75,.5,.25, 0)
#factor_propOfCR <- c(.02, .04, .06, .08, .10, .12, .14, .16, .18, .20,
#                     .22, .24, .26, .28, .30, .32, .34, .36, .38, .40)
# no of replications
#r <- 1:1000
# if(N != 1000) {print("Achtung!!!!!!!!!!!!!!!!! Geplante SP-Groesse ist 1000")}



# matrix with simulation conditions
design <- expand.grid("replication" = r,
                      "propOfCR" =  factor_propOfCR,
                      "randomOrder" = factor_randomOrder,
                      "stroRepOrder" = factor_strongRepOrder,
                      "modRepOrder" = factor_moderateRepOrder)
head(design)

design$relevant <- NA
# proportion can not be larger than 1
# no careless responding is also not reasonable
# either strong or moderate repetions of rank orders
design$relevant <- ifelse(
  rowSums(design[, c("randomOrder", "stroRepOrder", "modRepOrder")])!= 1  |
  rowSums(design[, c("randomOrder", "stroRepOrder", "modRepOrder")]) == 0 |
  (design[, "stroRepOrder"]>0 & design[,"modRepOrder"]>0),
                          0, 1)
design <- design[design$relevant == 1, ]
design$relevant <- NULL


# simulation seed
design$simSeed <- paste0("318", 1:nrow(design))
design$simSeed

####------------------ start simulation -------------------####


#test

no_ranks2sample <- 20
ranks2sample <- ro
m_out <- matrix(data = NA, nrow = 20, ncol = b*m)

m_CRss_resp <- matrix(data = NA, nrow = N*prop*design[d, 2], ncol = b*m)
m_CRms_resp <- matrix(data = NA, nrow = N*prop*design[d, 3], ncol = b*m)



matrixtest <- matrix(data = c("132", "123", "321", "111", "222", "333"),
                     2, 3, TRUE)
test <- apply(matrixtest, MARGIN = 1, FUN = ranks2string)
t(test)
###