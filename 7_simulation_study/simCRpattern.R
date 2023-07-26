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

# function: ranks2num ----
#'
#' @title ranks2num
#' 
#' @description casts strings to numeric output
#' @param rankOrder string with three digits (e.g., "123")
#' @return matrix with numeric output (e.g., 1 2 3) 
#' 
ranks2num<- function(rankOrder){
  out <- rbind(as.numeric(substr(rankOrder, 1, 1)),
               as.numeric(substr(rankOrder, 2, 2)),
               as.numeric(substr(rankOrder, 3, 3)))
  return(out)
}

# function: sampRepRank ----
#'
#' @title sampRepRank
#' 
#' @description with this function different types of careless responding 
#' pattern can be modeled. The type of careless responding pattern depends
#' on the parameter no_rep.
#' 
#' @param m_out output matrix
#' @param no_rep number of repetitions per rank (1 -> randomOrder, 5 -> moderateRepOrder, 
#' 20 -> strongRepOrder)
#' @param ro rank orders
#' @param ro_prob probabilities of the rank orders
#' @param b number of blocks in the questionnaire
#' 
#' @return matrix containing the different types of careless response pattern
#'

sampRepRank <- function(m_out, no_rep, ro, ro_prob, b){             

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

# function: se ----
#' 
#' 
#' 
#' 
#' 
se <- function(data, cr.index){
  tp <- nrow(data[data$CR==1 & data[, cr.index]==1, ])
  fn <- nrow(data[data$CR==1 & data[, cr.index]==0, ])
  
  return(tp/(tp+fn))
}

sp <- function(data, cr.index){
  tn <- nrow(data[data$CR==0 & data[, cr.index]==0, ])
  fp <- nrow(data[data$CR==0 & data[, cr.index]==1, ])
  
  return(tn/(tn+fp))
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

# careless responding indices
index <- c("cs", "md", "lom", "loa", "so", "tv")

# cut-off values
df_cf <- data.frame(cs = c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90),
                    lom = c(1, 2, 3, 4, 5, 6, 7),
                    loa = c(0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00),      
                    so = c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
                    tv = c(0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95))

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

unlink("7_simulation_study/TIRT/*", recursive=TRUE)
unlink("7_simulation_study/replications/*", recursive=TRUE)

res <- foreach(j = 1:nrow(design), 
               .inorder = FALSE,
               .packages = c("CRinMFC", "TirtAutomation", "MFCblockInfo"),
               .errorhandling = "pass")%dopar%{
                 
                 # select condition
                 con <- NULL
                 con <- design[j, ]
                 
                 # path
                 folder_rj <- paste0("7_simulation_study/TIRT/r-", j)
                 
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
                 dir.create(folder_rj)
                 
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
                   m_CR_ro <- sampRepRank(m_CR_ro, 
                                          no_rep = 1, 
                                          ro = ro, 
                                          ro_prob = ro_prob, 
                                          b = b)
                 }
                 
                 if(con[, "modRepOrder"] != 0){
                   
                   # size of sub-sample
                   n_modRepOrder <- N*(con[, "propOfCR"])*(con[, "modRepOrder"])
                   
                   # matrix for responses
                   m_CR_mro <- matrix(NA, nrow = n_modRepOrder, ncol = b*m)
                   
                   # sample rank orders
                   m_CR_mro <- sampRepRank(m_CR_mro, 
                                           no_rep = 5,
                                           ro = ro, 
                                           ro_prob = ro_prob, 
                                           b = b)
                 }
                 
                 if(con[, "stroRepOrder"] != 0){
                   
                   # size of sub-sample
                   n_stroRepOrder <- N*(con[, "propOfCR"])*(con[, "stroRepOrder"])
                   
                   # matrix for responses
                   m_CR_sro <- matrix(NA, nrow = n_stroRepOrder, ncol = b*m)
                   
                   # sample rank orders
                   m_CR_sro <- sampRepRank(m_CR_sro, 
                                           no_rep = 20,
                                           ro = ro, 
                                           ro_prob = ro_prob, 
                                           b = b)
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
                 multiplex::write.dat(bft, folder_rj)
                 
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
                     file=paste0(folder_rj, "/tirt-bft.inp"))
                 
                 # run TIRT
                 MplusAutomation::runModels(target = folder_rj)
                 
                 # read TIRT results an check for warnings
                 tirt <- MplusAutomation::readModels(target = folder_rj)
                 
                 # import parameters from mplus output
                 bt <- MplusAutomation::readModels(folder_rj, 
                                                   what=c("parameters", "sampstat"))
                 bt.pars <- bt$parameters$unstandardized
                 
                 # factor scores
                  m_theta.bt <- tryCatch(
                   expr = {
                     read.table(paste0(folder_rj, "/fs.dat"), header = FALSE, 
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
                    for(i in 1:nrow(df_cf)) {
                      dcri[,paste0("lom.cr", i)] <- ifelse(dcri$lom > df_cf[i, "lom"], 1, 0)
                      dcri[,paste0("loa.cr", i)] <- ifelse(dcri$loa > df_cf[i, "loa"], 1, 0)
                      dcri[,paste0("so.cr", i)] <- ifelse(dcri$so > df_cf[i, "so"], 1, 0)
                      dcri[,paste0("tv.cr", i)] <- ifelse(dcri$tv < df_cf[i, "tv"], 1, 0)
                    }
                    
                    
                    # compute performance measures according to Lalkhen & McCluskey (2008) ----
                    # se - sensitivity
                    # sp - specificity
                    
                    rep_result <- data.frame(seed = design[j, "simSeed"],
                                             rep = design[j, "replication"],
                                             cs_mean = mean(dcri$cs),
                                             md_mean = mean(dcri$md),
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
                                             
                                             cs_se1 = NA,
                                             cs_se2 = NA,
                                             cs_se3 = NA,
                                             cs_se4 = NA,
                                             cs_se5 = NA,
                                             cs_se6 = NA,
                                             cs_se7 = NA,
                                             
                                             cs_sp1 = NA,
                                             cs_sp2 = NA,
                                             cs_sp3 = NA,
                                             cs_sp4 = NA,
                                             cs_sp5 = NA,
                                             cs_sp6 = NA,
                                             cs_sp7 = NA,
                                             
                                             md_se = NA,
                                             md_sp = NA,
                                             
                                             lom_se1 = se(dcri, "lom.cr1"),
                                             lom_se2 = se(dcri, "lom.cr2"),
                                             lom_se3 = se(dcri, "lom.cr3"),
                                             lom_se4 = se(dcri, "lom.cr4"),
                                             lom_se5 = se(dcri, "lom.cr5"),
                                             lom_se6 = se(dcri, "lom.cr6"),
                                             lom_se7 = se(dcri, "lom.cr7"),
                                             
                                             lom_sp1 = sp(dcri, "lom.cr1"),
                                             lom_sp2 = sp(dcri, "lom.cr2"),
                                             lom_sp3 = sp(dcri, "lom.cr3"),
                                             lom_sp4 = sp(dcri, "lom.cr4"),
                                             lom_sp5 = sp(dcri, "lom.cr5"),
                                             lom_sp6 = sp(dcri, "lom.cr6"),
                                             lom_sp7 = sp(dcri, "lom.cr7"),
                                             
                                             loa_se1 = se(dcri, "loa.cr1"),
                                             loa_se2 = se(dcri, "loa.cr2"),
                                             loa_se3 = se(dcri, "loa.cr3"),
                                             loa_se4 = se(dcri, "loa.cr4"),
                                             loa_se5 = se(dcri, "loa.cr5"),
                                             loa_se6 = se(dcri, "loa.cr6"),
                                             loa_se7 = se(dcri, "loa.cr7"),
                                             
                                             loa_sp1 = sp(dcri, "loa.cr1"),
                                             loa_sp2 = sp(dcri, "loa.cr2"),
                                             loa_sp3 = sp(dcri, "loa.cr3"),
                                             loa_sp4 = sp(dcri, "loa.cr4"),
                                             loa_sp5 = sp(dcri, "loa.cr5"),
                                             loa_sp6 = sp(dcri, "loa.cr6"),
                                             loa_sp7 = sp(dcri, "loa.cr7"),
                                             
                                             so_se1 = se(dcri, "so.cr1"),
                                             so_se2 = se(dcri, "so.cr2"),
                                             so_se3 = se(dcri, "so.cr3"),
                                             so_se4 = se(dcri, "so.cr4"),
                                             so_se5 = se(dcri, "so.cr5"),
                                             so_se6 = se(dcri, "so.cr6"),
                                             so_se7 = se(dcri, "so.cr7"),
                                             
                                             so_sp1 = sp(dcri, "so.cr1"),
                                             so_sp2 = sp(dcri, "so.cr2"),
                                             so_sp3 = sp(dcri, "so.cr3"),
                                             so_sp4 = sp(dcri, "so.cr4"),
                                             so_sp5 = sp(dcri, "so.cr5"),
                                             so_sp6 = sp(dcri, "so.cr6"),
                                             so_sp7 = sp(dcri, "so.cr7"),
                                             
                                             tv_se1 = se(dcri, "tv.cr1"),
                                             tv_se2 = se(dcri, "tv.cr2"),
                                             tv_se3 = se(dcri, "tv.cr3"),
                                             tv_se4 = se(dcri, "tv.cr4"),
                                             tv_se5 = se(dcri, "tv.cr5"),
                                             tv_se6 = se(dcri, "tv.cr6"),
                                             tv_se7 = se(dcri, "tv.cr7"),
                                             
                                             tv_sp1 = sp(dcri, "tv.cr1"),
                                             tv_sp2 = sp(dcri, "tv.cr2"),
                                             tv_sp3 = sp(dcri, "tv.cr3"),
                                             tv_sp4 = sp(dcri, "tv.cr4"),
                                             tv_sp5 = sp(dcri, "tv.cr5"),
                                             tv_sp6 = sp(dcri, "tv.cr6"),
                                             tv_sp7 = sp(dcri, "tv.cr7"),
                                             error = "no_fs")
                    saveRDS(rep_result, file = paste0("7_simulation_study/replications/", 
                                                      "repRes", design[j, "simSeed"], ".RDS"))
                    
                    # delete folder with TIRT input and output data
                    unlink(folder_rj, recursive = TRUE)
                    
                    
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
                 dcri$md.cr <- ifelse(abs(dcri$md) > qchisq(p = 1-0.05, df = 60), 1, 0)
                 
                 for(i in 1:nrow(df_cf)) {
                   dcri[,paste0("cs.cr", i)] <- ifelse(dcri$cs < df_cf[i, "cs"], 1, 0)
                   dcri[,paste0("lom.cr", i)] <- ifelse(dcri$lom > df_cf[i, "lom"], 1, 0)
                   dcri[,paste0("loa.cr", i)] <- ifelse(dcri$loa > df_cf[i, "loa"], 1, 0)
                   dcri[,paste0("so.cr", i)] <- ifelse(dcri$so > df_cf[i, "so"], 1, 0)
                   dcri[,paste0("tv.cr", i)] <- ifelse(dcri$tv < df_cf[i, "tv"], 1, 0)
                 }
                 
                 
                 # compute performance measures according to Lalkhen & McCluskey (2008) ----
                 # se - sensitivity
                 # sp - specificity
                 
                 md.tp <- nrow(dcri[dcri$CR==1 & dcri$md.cr==1, ])
                 md.tn <- nrow(dcri[dcri$CR==0 & dcri$md.cr==0, ])
                 md.fp <- nrow(dcri[dcri$CR==0 & dcri$md.cr==1, ])
                 md.fn <- nrow(dcri[dcri$CR==1 & dcri$md.cr==0, ])
                 md.se <- md.tp/(md.tp+md.fn)
                 md.sp <- md.tn/(md.tn+md.fp)
                 
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
                                          
                                          cs_se1 = se(dcri, "cs.cr1"),
                                          cs_se2 = se(dcri, "cs.cr2"),
                                          cs_se3 = se(dcri, "cs.cr3"),
                                          cs_se4 = se(dcri, "cs.cr4"),
                                          cs_se5 = se(dcri, "cs.cr5"),
                                          cs_se6 = se(dcri, "cs.cr6"),
                                          cs_se7 = se(dcri, "cs.cr7"),
                                          
                                          cs_sp1 = sp(dcri, "cs.cr1"),
                                          cs_sp2 = sp(dcri, "cs.cr2"),
                                          cs_sp3 = sp(dcri, "cs.cr3"),
                                          cs_sp4 = sp(dcri, "cs.cr4"),
                                          cs_sp5 = sp(dcri, "cs.cr5"),
                                          cs_sp6 = sp(dcri, "cs.cr6"),
                                          cs_sp7 = sp(dcri, "cs.cr7"),
                                          
                                          md_se = md.se,
                                          md_sp = md.sp,
                                          
                                          lom_se1 = se(dcri, "lom.cr1"),
                                          lom_se2 = se(dcri, "lom.cr2"),
                                          lom_se3 = se(dcri, "lom.cr3"),
                                          lom_se4 = se(dcri, "lom.cr4"),
                                          lom_se5 = se(dcri, "lom.cr5"),
                                          lom_se6 = se(dcri, "lom.cr6"),
                                          lom_se7 = se(dcri, "lom.cr7"),
                                          
                                          lom_sp1 = sp(dcri, "lom.cr1"),
                                          lom_sp2 = sp(dcri, "lom.cr2"),
                                          lom_sp3 = sp(dcri, "lom.cr3"),
                                          lom_sp4 = sp(dcri, "lom.cr4"),
                                          lom_sp5 = sp(dcri, "lom.cr5"),
                                          lom_sp6 = sp(dcri, "lom.cr6"),
                                          lom_sp7 = sp(dcri, "lom.cr7"),
                                          
                                          loa_se1 = se(dcri, "loa.cr1"),
                                          loa_se2 = se(dcri, "loa.cr2"),
                                          loa_se3 = se(dcri, "loa.cr3"),
                                          loa_se4 = se(dcri, "loa.cr4"),
                                          loa_se5 = se(dcri, "loa.cr5"),
                                          loa_se6 = se(dcri, "loa.cr6"),
                                          loa_se7 = se(dcri, "loa.cr7"),
                                          
                                          loa_sp1 = sp(dcri, "loa.cr1"),
                                          loa_sp2 = sp(dcri, "loa.cr2"),
                                          loa_sp3 = sp(dcri, "loa.cr3"),
                                          loa_sp4 = sp(dcri, "loa.cr4"),
                                          loa_sp5 = sp(dcri, "loa.cr5"),
                                          loa_sp6 = sp(dcri, "loa.cr6"),
                                          loa_sp7 = sp(dcri, "loa.cr7"),
                                          
                                          so_se1 = se(dcri, "so.cr1"),
                                          so_se2 = se(dcri, "so.cr2"),
                                          so_se3 = se(dcri, "so.cr3"),
                                          so_se4 = se(dcri, "so.cr4"),
                                          so_se5 = se(dcri, "so.cr5"),
                                          so_se6 = se(dcri, "so.cr6"),
                                          so_se7 = se(dcri, "so.cr7"),
                                          
                                          so_sp1 = sp(dcri, "so.cr1"),
                                          so_sp2 = sp(dcri, "so.cr2"),
                                          so_sp3 = sp(dcri, "so.cr3"),
                                          so_sp4 = sp(dcri, "so.cr4"),
                                          so_sp5 = sp(dcri, "so.cr5"),
                                          so_sp6 = sp(dcri, "so.cr6"),
                                          so_sp7 = sp(dcri, "so.cr7"),
                                          
                                          tv_se1 = se(dcri, "tv.cr1"),
                                          tv_se2 = se(dcri, "tv.cr2"),
                                          tv_se3 = se(dcri, "tv.cr3"),
                                          tv_se4 = se(dcri, "tv.cr4"),
                                          tv_se5 = se(dcri, "tv.cr5"),
                                          tv_se6 = se(dcri, "tv.cr6"),
                                          tv_se7 = se(dcri, "tv.cr7"),
                                          
                                          tv_sp1 = sp(dcri, "tv.cr1"),
                                          tv_sp2 = sp(dcri, "tv.cr2"),
                                          tv_sp3 = sp(dcri, "tv.cr3"),
                                          tv_sp4 = sp(dcri, "tv.cr4"),
                                          tv_sp5 = sp(dcri, "tv.cr5"),
                                          tv_sp6 = sp(dcri, "tv.cr6"),
                                          tv_sp7 = sp(dcri, "tv.cr7"),
                                        
                                          error = NA)
                 saveRDS(rep_result, file = paste0("7_simulation_study/replications/", 
                                                   "repRes", design[j, "simSeed"], ".RDS"))
                 
                 # delete folder with TIRT input and output data
                 unlink(folder_rj, recursive = TRUE)
                 
               }

stopCluster(cl)

###

# combine all replications to one data frame and save the result
all_files <- list.files("7_simulation_study/replications", full.names = TRUE)
data <- lapply(all_files, readRDS)
dat <- do.call(rbind, data)

saveRDS(dat, "7_simulation_study/simRes230721.RDS")

###
