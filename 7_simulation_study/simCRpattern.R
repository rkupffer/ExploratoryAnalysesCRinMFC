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
                                           replace = TRUE,
                                           prob = ro_prob))
    
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
                                           replace = TRUE,
                                           prob = ro_prob))
    
    
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
                        replace = TRUE,
                        prob = ro_prob)
    
    # character to numeric ranks
    m_out_num <- ranks2num(rankOrders)
    
    # create output matrix and repeat the drwan rank order
    m_out <- matrix(data = cbind(replicate(no_rep, t(m_out_num))),
                    ncol = ncol(m_out), nrow = nrow(m_out),
                    byrow = FALSE)  
  }
  return(m_out)
}

#######################################################################
####------------------- simulation design -------------------------####

# fixed parameter ----

# sample and BFT questionnaire meta data ----
all_items <- read.csv2("2_All_Items_Coding.csv", 
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

# vector of rank order probabilities
ro_prob <- c(.25, rep(.15, 5))

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
int <- c(-1, 1)

# varying parameter ----

# design
p <- seq(0, 1, .25)
factor_randomOrder <- p
factor_strongRepOrder <- p
factor_moderateRepOrder <- p
factor_propOfCR <- seq(.02, .4, by = .05)

# replications
r <- 1:1000

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
# either strong or moderate repetitions of rank orders
design$relevant <- ifelse(
  rowSums(design[, c("randomOrder", "stroRepOrder", "modRepOrder")])!= 1  |
    rowSums(design[, c("randomOrder", "stroRepOrder", "modRepOrder")]) == 0 |
    (design[, "stroRepOrder"]>0 & design[,"modRepOrder"]>0),
  0, 1)
design <- design[design$relevant == 1, ]
design$relevant <- NULL

# simulation seed
design$simSeed <- paste0("318", 1:nrow(design))


###############################################################
###############################################################
####------------------ start simulation -------------------####

cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

unlink("TIRT/*", recursive=TRUE)
unlink("replications/*", recursive=TRUE)

res <- foreach(j = 1:nrow(design), 
               .inorder = FALSE,
               .packages = c("CRinMFC", "TirtAutomation", "MFCblockInfo",),
               .errorhandling = "pass")%dopar%{
                 
                 # select condition
                 con <- NULL
                 con <- design[j, ]
                 
                 # reset matrices 
                 m_Care_resp <- NULL
                 m_CR_ro <- NULL
                 m_CR_mro <- NULL
                 m_CR_sro <- NULL
                 traits <- NULL
                 rep_result <- NULL
                 
                 #set simulation seed
                 set.seed(design[j, "simSeed"])
                
                 
                 # simulation of careful/thoughtful responses ----
                 # sample size of careful subsample
                 n_care <- round(N*(1-con[, "propOfCR"]))
                 
                 # draw traits of n_care participants ----
                 traits <- mvtnorm::rmvnorm(n_care, v_mu, m_phi, method = "chol")
                 
                 # if no error occurs, create a folder for Mplus inp, out & data
                 dir.create(sprintf("TIRT/r-%i", j))
                 
                 # simulate item parameter ----
                 bft_items <- sim.items(m_dload, b, m, load, int)
                 
                 # simulation of the responses as ranks
                 resp <- sim.responses(traits, bft_items, m_dload, b, m, return.index = FALSE)
                 
                 # save as matrix
                 m_Care_resp <- resp$ranks
                 
                 
                 # simulation of different careless responding pattern ----
                 
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
                 multiplex::write.dat(bft, sprintf("TIRT/r-%i", j))
                 
                 # rename traits an variables in design load matrix
                 colnames(m_dload) <- c("TNeu", "TExt", "TOpe", "TAgr", "TCon")
                 rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("T","0",1:nrow(m_dload)), paste0("T",1:nrow(m_dload)))
                 
                 tirt_bft <- tirtMplusSyntax(design.load = m_dload, 
                                             names.pairs = NULL,
                                             item.short = "T",
                                             id.var = "ID",
                                             file.data = "bft.dat",
                                             title = paste0("bft", design[j, "simSeed"]), 
                                             out.command = "sampstat standardized;",
                                             fscores.file = "fs.dat",
                                             missings.as = "-99")
                 # save as input-file
                 cat(paste(tirt_bft, collapse="\n\n"), 
                     file=paste0(sprintf("TIRT/r-%i/tirt-bft.inp", j)))
                 
                 # run TIRT
                 MplusAutomation::runModels(target = sprintf("TIRT/r-%i", j))
                 
                 # read TIRT results an check for warnings
                 tirt <- MplusAutomation::readModels(target = sprintf("TIRT/r-%i", j))
                 
                 # import parameters from mplus output
                 bt <- MplusAutomation::readModels(sprintf("TIRT/r-%i", j), 
                                                   what=c("parameters", "sampstat"))
                 bt.pars <- bt$parameters$unstandardized
                 
                 # factor scores
                  m_theta.bt <- tryCatch(
                   expr = {
                     read.table(sprintf("TIRT/r-%i/fs.dat", j), header = FALSE, 
                                na.strings = "*")[,c(71,61,63,65,67,69)]
                   },
                   error = function(e){
                     print(e)
                     FALSE
                   },
                   warning = function(w){
                     print(w)
                     FALSE
                   }
                  )
                  
                  if(is.logical(m_theta.bt)){
                    # if no factor scores were computed, compute indices for
                    # which no factor scores are needed 
                    
                    # data frame of all CR indices ----
                    dcri <- data.frame(ID = bft[,1], 
                                       #was careless responding simulated?
                                       CR = c(rep(0, n_care), rep(1, nrow(m_CR_resp))))
                    
                    # LongOrderMax ----
                    blocks <- matrix(seq(1:I), m, b)
                    dat.b <- apply(blocks, 2, function(bn, d.r) 
                      apply(d.r[,bn], 1, paste, collapse=""), m_all)
                    
                    # calculate longOrder for each participant
                    m_lo <- longOrder(IDs = bft[, "ID"],
                                      d.b.o = dat.b,
                                      no.b = b)
                    
                    dcri$lom <- longOrderMax(m_lo)
                    
                    # LongOrderAvg ----
                    dcri$loa <- longOrderAvg(m_lo)
                    
                    # SameOrder ----
                    dcri$so <- sameOrder(dat.b, b)
                    
                    # Triplet Variance ----
                    dcri$tv <- tripletVariance(dat.b)
                    
                    # apply cut-off values
                    dcri$lom.cr <- ifelse(dcri$lom > 7, 1, 0)
                    dcri$loa.cr <- ifelse(dcri$loa > 2, 1, 0)
                    dcri$so.cr <- ifelse(dcri$so > .4, 1, 0)
                    dcri$tv.cr <- ifelse(dcri$tv < .7, 1, 0)
                    
                    # compute performance measures according to Lalkhen & McCluskey (2008) ----
                    # se - sensitivity
                    # sp - specificity
                    # ppv - positive predictive value
                    # npv - negative predictive value
                    
                    lom.tp <- nrow(dcri[dcri$CR==1 & dcri$lom.cr==1, ])
                    lom.tn <- nrow(dcri[dcri$CR==0 & dcri$lom.cr==0, ])
                    lom.fp <- nrow(dcri[dcri$CR==0 & dcri$lom.cr==1, ])
                    lom.fn <- nrow(dcri[dcri$CR==1 & dcri$lom.cr==0, ])
                    lom.se <- lom.tp/(lom.tp+lom.fn)
                    lom.sp <- lom.tn/(lom.tn+lom.fp)
                    lom.ppv <- lom.tp/(lom.tp+lom.fp)
                    lom.npv <- lom.tn/(lom.tn+lom.fn)
                    
                    loa.tp <- nrow(dcri[dcri$CR==1 & dcri$loa.cr==1, ])
                    loa.tn <- nrow(dcri[dcri$CR==0 & dcri$loa.cr==0, ])
                    loa.fp <- nrow(dcri[dcri$CR==0 & dcri$loa.cr==1, ])
                    loa.fn <- nrow(dcri[dcri$CR==1 & dcri$loa.cr==0, ])
                    loa.se <- loa.tp/(loa.tp+loa.fn)
                    loa.sp <- loa.tn/(loa.tn+loa.fp)
                    loa.ppv <- loa.tp/(loa.tp+loa.fp)
                    loa.npv <- loa.tn/(loa.tn+loa.fn)
                    
                    so.tp <- nrow(dcri[dcri$CR==1 & dcri$so.cr==1, ])
                    so.tn <- nrow(dcri[dcri$CR==0 & dcri$so.cr==0, ])
                    so.fp <- nrow(dcri[dcri$CR==0 & dcri$so.cr==1, ])
                    so.fn <- nrow(dcri[dcri$CR==1 & dcri$so.cr==0, ])
                    so.se <- so.tp/(so.tp+so.fn)
                    so.sp <- so.tn/(so.tn+so.fp)
                    so.ppv <- so.tp/(so.tp+so.fp)
                    so.npv <- so.tn/(so.tn+so.fn)
                    
                    tv.tp <- nrow(dcri[dcri$CR==1 & dcri$tv.cr==1, ])
                    tv.tn <- nrow(dcri[dcri$CR==0 & dcri$tv.cr==0, ])
                    tv.fp <- nrow(dcri[dcri$CR==0 & dcri$tv.cr==1, ])
                    tv.fn <- nrow(dcri[dcri$CR==1 & dcri$tv.cr==0, ])
                    tv.se <- tv.tp/(tv.tp+tv.fn)
                    tv.sp <- tv.tn/(tv.tn+tv.fp)
                    tv.ppv <- tv.tp/(tv.tp+tv.fp)
                    tv.npv <- tv.tn/(tv.tn+tv.fn)
                    
                    rep_result <- data.frame(seed = design[j, "simSeed"],
                                             rep = design[j, "replication"],
                                             cs_mean = NA,
                                             md_mean = NA,
                                             lom_mean = mean(dcri$lom),
                                             loa_mean = mean(dcri$loa),
                                             so_mean = mean(dcri$so),
                                             tv_mean = mean(dcri$tv),
                                             cs_sd = NA,
                                             md_sd = NA,
                                             lom_sd = sd(dcri$lom),
                                             loa_sd = sd(dcri$loa),
                                             so_sd = sd(dcri$so),
                                             tv_sd = sd(dcri$tv),
                                             cs_se = NA,
                                             cs_sp = NA,
                                             cs_ppv = NA,
                                             cs_npv = NA,
                                             md_se = NA,
                                             md_sp = NA,
                                             md_ppv = NA,
                                             md_npv = NA,
                                             lom_se = lom.se,
                                             lom_sp = lom.sp,
                                             lom_ppv = lom.ppv,
                                             lom_npv = lom.npv,
                                             loa_se = loa.se,
                                             loa_sp = loa.sp,
                                             loa_ppv = loa.ppv,
                                             loa_npv = loa.npv,
                                             so_se = so.se,
                                             so_sp = so.sp,
                                             so_ppv = so.ppv,
                                             so_npv = so.npv,
                                             tv_se = tv.se,
                                             tv_sp = tv.sp,
                                             tv_ppv = tv.ppv,
                                             tv_npv = tv.npv,
                                             error = "no_fs")
                    saveRDS(rep_result, file = paste0("replications/", 
                                                      "repRes", design[j, "simSeed"], ".RDS"))
                    
                    # delete folder with TIRT input and output data
                    unlink(paste0("TIRT/r-",j), recursive = TRUE)
                    
                    
                    return(paste0("no factor scores found ("
                                  , design[j, "simSeed"]), ")")
                  }
                  
                 
                 colnames(m_theta.bt)<- c("ID",paste0(c("N","E","O","A","C"),"_bt"))
                 
                 # data frame of all CR indices ----
                 dcri <- data.frame(ID = bft[,1], 
                                    #was careless responding simulated?
                                    CR = c(rep(0, n_care), rep(1, nrow(m_CR_resp))))
                 
                 # Consistency Score ----
                 dcri$cs <- consisScore(
                   quest.pars = bt.pars,
                   no.b = b,
                   no.traits = 5,
                   all.items = all_items,
                   quest = "BT",
                   m.theta = m_theta.bt[, 2:6],
                   d.bi = bft[, 2:61])
                 
                 # Mahalanobis Distance ----
                 cor.bt <- cor2mat(bt)
                 thresh.bt <- thresh2vec(bt)
                 d_center.bt <- sweep(bft[, 2:61], 2L, thresh.bt)
                 dcri$md <- mahalanobis(d_center.bt, center = FALSE, cov = cor.bt)
                 
                 # LongOrderMax ----
                 blocks <- matrix(seq(1:I), m, b)
                 dat.b <- apply(blocks, 2, function(bn, d.r) 
                   apply(d.r[,bn], 1, paste, collapse=""), m_all)
                 
                 # calculate longOrder for each participant
                 m_lo <- longOrder(IDs = bft[, "ID"],
                                   d.b.o = dat.b,
                                   no.b = b)
                 
                 dcri$lom <- longOrderMax(m_lo)
                 
                 # LongOrderAvg ----
                 dcri$loa <- longOrderAvg(m_lo)
                 
                 # SameOrder ----
                 dcri$so <- sameOrder(dat.b, b)
                 
                 # Triplet Variance ----
                 dcri$tv <- tripletVariance(dat.b)
                 
                 # apply cut-off values
                 dcri$cs.cr <- ifelse(dcri$cs < .64, 1, 0)
                 dcri$md.cr <- ifelse(abs(dcri$md) > qchisq(p = 1-0.05, df = 60), 1, 0)
                 dcri$lom.cr <- ifelse(dcri$lom > 7, 1, 0)
                 dcri$loa.cr <- ifelse(dcri$loa > 2, 1, 0)
                 dcri$so.cr <- ifelse(dcri$so > .4, 1, 0)
                 dcri$tv.cr <- ifelse(dcri$tv < .7, 1, 0)
                 
                 # compute performance measures according to Lalkhen & McCluskey (2008) ----
                 # se - sensitivity
                 # sp - specificity
                 # ppv - positive predictive value
                 # npv - negative predictive value
                 
                 # prefix (e.g., "cs.") is added to avoid errors in computation in case an index
                 # could not be computed
                 cs.tp <- nrow(dcri[dcri$CR==1 & dcri$cs.cr==1, ])
                 cs.tn <- nrow(dcri[dcri$CR==0 & dcri$cs.cr==0, ])
                 cs.fp <- nrow(dcri[dcri$CR==0 & dcri$cs.cr==1, ])
                 cs.fn <- nrow(dcri[dcri$CR==1 & dcri$cs.cr==0, ])
                 cs.se <- cs.tp/(cs.tp+cs.fn)
                 cs.sp <- cs.tn/(cs.tn+cs.fp)
                 cs.ppv <- cs.tp/(cs.tp+cs.fp)
                 cs.npv <- cs.tn/(cs.tn+cs.fn)
                 
                 md.tp <- nrow(dcri[dcri$CR==1 & dcri$md.cr==1, ])
                 md.tn <- nrow(dcri[dcri$CR==0 & dcri$md.cr==0, ])
                 md.fp <- nrow(dcri[dcri$CR==0 & dcri$md.cr==1, ])
                 md.fn <- nrow(dcri[dcri$CR==1 & dcri$md.cr==0, ])
                 md.se <- md.tp/(md.tp+md.fn)
                 md.sp <- md.tn/(md.tn+md.fp)
                 md.ppv <- md.tp/(md.tp+md.fp)
                 md.npv <- md.tn/(md.tn+md.fn)
                 
                 lom.tp <- nrow(dcri[dcri$CR==1 & dcri$lom.cr==1, ])
                 lom.tn <- nrow(dcri[dcri$CR==0 & dcri$lom.cr==0, ])
                 lom.fp <- nrow(dcri[dcri$CR==0 & dcri$lom.cr==1, ])
                 lom.fn <- nrow(dcri[dcri$CR==1 & dcri$lom.cr==0, ])
                 lom.se <- lom.tp/(lom.tp+lom.fn)
                 lom.sp <- lom.tn/(lom.tn+lom.fp)
                 lom.ppv <- lom.tp/(lom.tp+lom.fp)
                 lom.npv <- lom.tn/(lom.tn+lom.fn)
                 
                 loa.tp <- nrow(dcri[dcri$CR==1 & dcri$loa.cr==1, ])
                 loa.tn <- nrow(dcri[dcri$CR==0 & dcri$loa.cr==0, ])
                 loa.fp <- nrow(dcri[dcri$CR==0 & dcri$loa.cr==1, ])
                 loa.fn <- nrow(dcri[dcri$CR==1 & dcri$loa.cr==0, ])
                 loa.se <- loa.tp/(loa.tp+loa.fn)
                 loa.sp <- loa.tn/(loa.tn+loa.fp)
                 loa.ppv <- loa.tp/(loa.tp+loa.fp)
                 loa.npv <- loa.tn/(loa.tn+loa.fn)
                 
                 so.tp <- nrow(dcri[dcri$CR==1 & dcri$so.cr==1, ])
                 so.tn <- nrow(dcri[dcri$CR==0 & dcri$so.cr==0, ])
                 so.fp <- nrow(dcri[dcri$CR==0 & dcri$so.cr==1, ])
                 so.fn <- nrow(dcri[dcri$CR==1 & dcri$so.cr==0, ])
                 so.se <- so.tp/(so.tp+so.fn)
                 so.sp <- so.tn/(so.tn+so.fp)
                 so.ppv <- so.tp/(so.tp+so.fp)
                 so.npv <- so.tn/(so.tn+so.fn)
                 
                 tv.tp <- nrow(dcri[dcri$CR==1 & dcri$tv.cr==1, ])
                 tv.tn <- nrow(dcri[dcri$CR==0 & dcri$tv.cr==0, ])
                 tv.fp <- nrow(dcri[dcri$CR==0 & dcri$tv.cr==1, ])
                 tv.fn <- nrow(dcri[dcri$CR==1 & dcri$tv.cr==0, ])
                 tv.se <- tv.tp/(tv.tp+tv.fn)
                 tv.sp <- tv.tn/(tv.tn+tv.fp)
                 tv.ppv <- tv.tp/(tv.tp+tv.fp)
                 tv.npv <- tv.tn/(tv.tn+tv.fn)
                 
                 rep_result <- data.frame(seed = design[j, "simSeed"],
                                          rep = design[j, "replication"],
                                          cs_mean = mean(dcri$cs),
                                          md_mean = mean(dcri$md),
                                          lom_mean = mean(dcri$lom),
                                          loa_mean = mean(dcri$loa),
                                          so_mean = mean(dcri$so),
                                          tv_mean = mean(dcri$tv),
                                          cs_sd = sd(dcri$cs),
                                          md_sd = sd(dcri$md),
                                          lom_sd = sd(dcri$lom),
                                          loa_sd = sd(dcri$loa),
                                          so_sd = sd(dcri$so),
                                          tv_sd = sd(dcri$tv),
                                          cs_se = cs.se,
                                          cs_sp = cs.sp,
                                          cs_ppv = cs.ppv,
                                          cs_npv = cs.npv,
                                          md_se = md.se,
                                          md_sp = md.sp,
                                          md_ppv = md.ppv,
                                          md_npv = md.npv,
                                          lom_se = lom.se,
                                          lom_sp = lom.sp,
                                          lom_ppv = lom.ppv,
                                          lom_npv = lom.npv,
                                          loa_se = loa.se,
                                          loa_sp = loa.sp,
                                          loa_ppv = loa.ppv,
                                          loa_npv = loa.npv,
                                          so_se = so.se,
                                          so_sp = so.sp,
                                          so_ppv = so.ppv,
                                          so_npv = so.npv,
                                          tv_se = tv.se,
                                          tv_sp = tv.sp,
                                          tv_ppv = tv.ppv,
                                          tv_npv = tv.npv,
                                          error = NA)
                 saveRDS(rep_result, file = paste0("replications/", 
                                                   "repRes", design[j, "simSeed"], ".RDS"))
                 
                 # delete folder with TIRT input and output data
                 unlink(paste0("TIRT/r-",j), recursive = TRUE)
                 
               }

stopCluster(cl)

###

# combine all replications to one data frame and save the result

data <- lapply(list.files("replications", full.names = TRUE),
               readRDS)
dat <- do.call(rbind, data)

saveRDS(dat, "simRes7200.RDS")

###
