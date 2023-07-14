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

sampRepRank <- function(m_out,           # output matrix
                        no_rep){         # number of repetitions per rank
                                         # 1 - randomOrder
                                         # 5 - moderateRepOrder
                                         # 20- strongRepOrder
  
  # no of ranks drawn (per participant)
  no_ranks2sample <- b/no_rep
  
  # random selection of rank orders - randomOrder  ------------------------
  if(no_rep == 1){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(rankOrders, 
                        MARGIN = 2, 
                        function(x) sample(ro, 
                                           size = nrow(m_out), 
                                           replace = TRUE))
    
    # character to numeric ranks
    m_out_num <- matrix(data = NA, nrow = nrow(m_out), ncol = no_ranks2sample*3)
    
    m_out_num <- apply(rankOrders,
                       MARGIN = 1,
                       FUN = ranks2num)
    
    m_out <- t(m_out_num)
    
  }
  
  # moderate repetition of rank orders - modRepOrder ------------------------
  if(no_rep == 5){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(rankOrders, 
                        MARGIN = 1:2, 
                        function(x) sample(ro, 
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
  
  # strong repetition of rank orders - stroRepOrder  ------------------------
  if(no_rep == 20){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(ro, 
                        MARGIN = 2, 
                        FUN = sample,
                        size = nrow(m_out), 
                        replace = TRUE)
    
    # character to numeric ranks
    m_out_num <- ranks2num(rankOrders)
    
    # create output matrix and repeat the drwan rank order
    m_out <- matrix(data = cbind(replicate(no_rep, t(m_out_num))),
                    ncol = ncol(m_out), nrow = nrow(m_out),
                    byrow = FALSE)  
  }
  return(m_out)
}

####------------------- simulation design -------------------------####

# fixed parameter ----

# sample and BFT questionnaire meta data ----
all_items <- read.csv2("../ExploratoryAnalysesCRinMFC/2_All_Items_Coding.csv", 
                       sep = ",")[,2:8]

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "BT", no.traits = 5)

# number of blocks in the questionnaire
b <- 20

# number of items per block
m <- 3

# total number of items
I <- b*m

# design matrix A: rows - items, cols - pairwise comparisons 
A <- TirtAutomation::designMatrix(no.b = b)

# sample size
N <- 1000

# rank orders stored in a matrix
ro <- matrix(c("123", "132", "213", "231", "312", "321"), 6, byrow = TRUE)

# variance-covariance matrix (Anglim et al., 2020, p. 60)
m_phi <- matrix(data = c(1, -.49, -.19, -.17, -.48, 
                         -.49, 1, .30, .08, .19, 
                         -.19, .30, 1, .20, .14,
                         -.17, .08, .20, 1, .32,
                         -.48, .19, .14, .32, 1), 
                nrow = 5, ncol = 5, byrow = TRUE)

# vector with means fixed to 0
v_mu <- rep(0, 5)

# range of factor loadings
load <- c(.65, .8)

# range of intercepts
int <- c(-1, .1)

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
N <- 1000




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

i <- 101

####------------------ start simulation -------------------####

# cl <- makeCluster(10)
# registerDoParallel(cl)

for(i in 1:nrow(design)){
  
  #set simulation seed
  set.seed(design[i, "simSeed"])
  
  # select condition
  con <- NULL
  con <- design[i, ]
  
  # reset matrices 
  m_Care_resp <- NULL
  m_CR_ro <- NULL
  m_CR_mro <- NULL
  m_CR_sro <- NULL
  
  # delete content of the previous TIRT folder
  unlink("7_simulation_study/TIRT/*")
  
  # simulation of careful/thoughtful responses -------------------------------
  # sample size of careful subsample
  n_care <- N*(1-con[, "propOfCR"])
  
  # simulate item parameter ----
  bft_items <- sim.items(m_dload, b, m, load, int)
  
  # draw traits of n_care participants
  traits <- mvtnorm::rmvnorm(n_care, v_mu, m_phi, method = "chol")
  
  # simulation of the responses as ranks
  resp <- sim.responses(traits, bft_items, m_dload, b, m, return.index = FALSE)
  
  # save as matrix
  m_Care_resp <- resp$ranks
  
  # simulation of different careless responding pattern ----------------------
  
  if(con[, "randomOrder"] != 0){
    
    # size of sub-sample
    n_randomOrder <- N*(con[, "propOfCR"])*(con[, "randomOrder"])
    
    # matrix for responses
    m_CR_ro <- matrix(NA, nrow = n_randomOrder, ncol = b*m)
    
    # sample rank orders
    m_CR_ro <- sampRepRank(m_CR_ro, no_rep = 1)
  }
  
  if(con[, "modRepOrder"] != 0){
    
    # size of sub-sample
    n_modRepOrder <- N*(con[, "propOfCR"])*(con[, "modRepOrder"])
    
    # matrix for responses
    m_CR_mro <- matrix(NA, nrow = n_modRepOrder, ncol = b*m)
    
    # sample rank orders
    m_CR_mro <- sampRepRank(m_CR_mro, no_rep = 5)
  }
  
  if(con[, "stroRepOrder"] != 0){
    
    # size of sub-sample
    n_stroRepOrder <- N*(con[, "propOfCR"])*(con[, "stroRepOrder"])
    
    # matrix for responses
    m_CR_sro <- matrix(NA, nrow = n_stroRepOrder, ncol = b*m)
    
    # sample rank orders
    m_CR_sro <- sampRepRank(m_CR_sro, no_rep = 20)
  }
  
  # combine all response pattern in this condition to one matrix -------------
  m_CR_resp <- NULL
  
  if(!is.null(m_CR_ro) && !is.null(m_CR_mro)){
    m_CR_resp <- rbind(m_CR_ro, m_CR_mro)
  } else if(!is.null(m_CR_ro) && !is.null(m_CR_sro)){
    m_CR_resp <- rbind(m_CR_ro, m_CR_sro)
  } else if(!is.null(m_CR_ro) && is.null(m_CR_mro) && is.null(m_CR_sro)){
    m_CR_resp <- m_CR_ro
  } else if(is.null(m_CR_ro) && !is.null(m_CR_mro)){
    m_CR_resp <- m_CR_mro
  } else{
    m_CR_resp <- m_CR_sro
  }
  m_all <- NULL
  m_all <- rbind(m_Care_resp, m_CR_resp)
  
  # apply MFC careless responding indices ------------------------------------
  ID <- 1:nrow(m_all)
  
  # recode rank order to binary outcomes
  bft <- reTri(as.data.frame(m_all), "b")
  bft <- cbind(ID, bft)
  
  # save binary outcomes
  multiplex::write.dat(bft, "7_simulation_study/TIRT")
  
  # rename traits an variables in design load matrix
  colnames(m_dload) <- c("TNeu", "TExt", "TOpe", "TAgr", "TCon")
  rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("T","0",1:nrow(m_dload)), paste0("T",1:nrow(m_dload)))
  
  tirt_bft <- tirtMplusSyntax(design.load = m_dload, 
                              names.pairs = NULL,
                              item.short = "T",
                              id.var = "ID",
                              file.data = "bft.dat", 
                              title = paste0("bft", design[i, "simSeed"]), 
                              out.command = "sampstat standardized;",
                              fscores.file = paste0("fs.dat"),
                              missings.as = "-99")
  # save as input-file
  cat(paste(tirt_bft, collapse="\n\n"), 
      file=paste0("7_simulation_study/TIRT/tirt-bft.inp"))
  
  # run TIRT
  MplusAutomation::runModels(target = "7_simulation_study/TIRT")
  
  # read TIRT results an check for warnings
  tirt <- MplusAutomation::readModels(target = "7_simulation_study/TIRT")
  tirt$warnings
  tirt$errors
  
  # -------------- Consistency Score -----------------------------------
  # import parameters from mplus output
  bt <- MplusAutomation::readModels("7_simulation_study/TIRT", what=c("parameters", "sampstat"))
  bt.pars <- bt$parameters$unstandardized
  
  # factor scores
  m_theta.bt <- read.table("7_simulation_study/TIRT/fs.dat", header = FALSE, na.strings = "*")[,c(71,61,63,65,67,69)]
  colnames(m_theta.bt)<- c("ID",paste0(c("N","E","O","A","C"),"_bt"))
  
  # compute consistency score
  cs <- consisScore(
    quest.pars = bt.pars,
    no.b = b,
    no.traits = 5,
    all.items = all_items,
    quest = "BT",
    m.theta = m_theta.bt[, 2:6],
    d.bi = bft[, 2:61])
  
  # -------------- Mahalanobis Distance --------------------------------
  cor.bt <- cor2mat(bt)
  thresh.bt <- thresh2vec(bt)
  d_center.bt <- sweep(bft[, 2:61], 2L, thresh.bt)
  md <- mahalanobis(d_center.bt, center = FALSE, cov = cor.bt)
  
  # -------------- LongOrderMax ----------------------------------------
  blocks <- matrix(I, m, b)
  dat.b <- apply(blocks, 2, function(bn, d.r) 
    apply(d.r[,bn], 1, paste, collapse=""), m_all[, 2:61])
  
  # calculate longOrder for each participant
  m_lo <- longOrder(IDs = m_all[, "ID"],
                    d.b.o = dat.b,
                    no.b = b)
  
  lom <- longOrderMax(m_lo)
  
  # -------------- LongOrderAvg ----------------------------------------
  loa <- longOrderAvg(m_lo)
  
  # -------------- SameOrder -------------------------------------------
  so <- sameOrder(m_all[, 2:61], no.b = b)
  
  # -------------- Triplet Variance ------------------------------------
  tv <- tripletVariance(m_all[, 2:61])
  
  # --------------  table of all CR indices ----------------------------
  
  
} ### end of sim




###
