#### simulation and detection of different careless responding pattern  ####

# In this script different careless and careful responding pattern 
# will be simulated with varying proportions of careless respondents.

# packages ----
library(MFCblockInfo)
library(TirtAutomation)
library(CRinMFC)

rm(list = ls())


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

# function: repeat selected rank order ----

repeatRank <- function(rank, no_rep){
  e1 <- as.numeric(substr(rank, 1, 1))
  e2 <- as.numeric(substr(rank, 2, 2))
  e3 <- as.numeric(substr(rank, 3, 3))
  
  out <- rep(c(e1, e2, e3), no_rep)
  out
}

# function: sample and repeat rank orders ----

sampRepRank <- function(ranks2sample,    # list with rank orders as strings
                        m_out,           # output matrix
                        b,               # no of blocks
                        no_ranks2sample, # how many ranks will be drawn
                                         # 1 - ss, 4 - ms, 20 - rs
                        no_rep){         # number of repetitions per rank
  
  # strong swiping
  if(no_ranks2sample == 1){
    
    # draw ranks and store them in a list 
    rankOrders <- matrix(NA, nrow = nrow(m_out), ncol = no_ranks2sample)
    rankOrders <- apply(ranks2sample, MARGIN = 2, FUN = sample,
                        size = nrow(m_out), replace = TRUE)
    
    # character to numeric ranks
    e1 <- as.numeric(substr(rankOrders, 1, 1))
    e2 <- as.numeric(substr(rankOrders, 2, 2))
    e3 <- as.numeric(substr(rankOrders, 3, 3))
    
    # repeat rank order no_rep times
    m_out <- rep(cbind(e1, e2, e3), b)
  }
  
  # moderate swiping
  if(no_ranks2sample == b/5){
    
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
    
    m_out <- matrix(data = cbind(replicate(5, m_ro1),
                                 replicate(5, m_ro2),
                                 replicate(5, m_ro3),
                                 replicate(5, m_ro4)),
                    ncol = ncol(m_out), nrow = nrow(m_out),
                    byrow = FALSE)  
    
  }
}

# simulation design ----
'design <- data.frame(random = c(rep(c(0,.25,.5,.75,1),2)),
                     strongS = c(1,.75,.5,.25, rep(0, 6)),
                     moderateS = c(rep(0, 5), 1,.75,.5,.25, 0)
                     )'

design <- data.frame(random = c(rep(c(0,.5,1),2)),
                     strongS = c(1,.5, rep(0, 4)),
                     moderateS = c(rep(0, 3), 1, .5, 0)
)

'propOfCR <- c(.02, .04, .06, .08, .10, .12, .14, .16, .18, .20,
              .22, .24, .26, .28, .30, .32, .34, .36, .38, .40)'

propOfCR <- c(.02, .07, .12, .17, .22, .27, .32, .37)

# number of replications
r <- 1:10 # later 1000

# simulation seed
simSeed <- paste0("318", 1:(length(r)*nrow(design)))
simSeed

# start ----

#TODO: save response matrices with corresponding seeds and the info what kind
#of response pattern was modeled

#TODO: check the loop. Something is off with the sequence and seed iteration.

#TODO: rewrite code to functions.


for(i in 1:length(simSeed)){
  
  set.seed(simSeed[i])
  replSeed <- simSeed[i]
  
  for(p in 1:length(propOfCR)){
    
    prop <- NULL
    prop <- propOfCR[p]
    
    for(d in 1:nrow(design)){
      
      # select condition
      con <- NULL
      con <- design[d, ]
      
      
      # simulation of careful/thoughtful responses ----
      
      # item parameter simulation
      bft_items <- sim.items(m_dload, b, m, load, int)
      
      # draw traits of N persons
      traits <- mvtnorm::rmvnorm(N*(1-prop), v_mu, m_sigma)
      
      # simulation of the responses as ranks
      resp <- sim.responses(traits, bft_items, m_dload, b, m, 
                            return.index = FALSE)
      m_Care_resp <- resp$ranks
      
      # random responding
      if(design[d, 1]!=0){
        m_CRrr_resp <- matrix(data = NA, nrow = N*prop*design[d, 1], ncol = b*m)
        
        for(i in 1:nrow(m_CRrr_resp)){
          for(r in seq(1, to = b*m, by = m)){
            rank <- sample(ro, 1)
            m_CRrr_resp[i, c(r, r+1, r+2)] <- repeatRank(rank, no_rep = 1)
          }
        }
      }
      
      # strong swiping
      if(design[d, 2]!=0){
        m_CRss_resp <- matrix(data = NA, nrow = N*prop*design[d, 2], ncol = b*m)
        
        for(i in 1:nrow(m_CRss_resp)){
          rank <- sample(ro, 1)
          m_CRss_resp[i, ] <- repeatRank(rank, no_rep = b)
        }
        m_CRss_resp
      }
      
      # moderate swiping
      triplet_rep <- 5
      
      if(design[d, 3]!=0){
        m_CRms_resp <- matrix(data = NA, nrow = N*prop*design[d, 3], ncol = b*m)
        
        for(i in 1:nrow(m_CRms_resp)){
          for(r in 1:triplet_rep){
            rank <- sample(ro, 1)
            m_CRms_resp[i, 1:(triplet_rep*3)] <- repeatRank(rank, no_rep = triplet_rep)
          }
          for(r in triplet_rep+1:triplet_rep*2){
            rank <- sample(ro, 1)
            m_CRms_resp[i, 16:30] <- repeatRank(rank, no_rep = triplet_rep)
          }
          for(r in triplet_rep*2+1:triplet_rep*3){
            rank <- sample(ro, 1)
            m_CRms_resp[i, 31:45] <- repeatRank(rank, no_rep = triplet_rep)
          }
          for(r in triplet_rep*3+1:triplet_rep*4){
            rank <- sample(ro, 1)
            m_CRms_resp[i, 46:60] <- repeatRank(rank, no_rep = triplet_rep)
          }
        }
      }
      
      
      
    }
  }
  
  #save as text doc with sim seed
  
}









###