######### additional functions ############
# ------------------------------- round2() -----------------------------------#
# round .5 up 

# source: http://andrewlandgraf.com/2012/06/15/rounding-in-r/
# author: Andrew J. Landgraf
#    and also discussed in an answer to a question on stackoverflow:
#    https://stackoverflow.com/questions/12688717/round-up-from-5-in-r

round2 = function(x, n) {
  
  posneg <- sign(x) #store sign
  z <- abs(x)*10^n  #move comma n decimals to the right
  z <- z + 0.5      #add 5 to (n+1)-th decimal
  z <- trunc(z)     #cut of other decimals
  z <- z/10^n       #move comma back
  z*posneg          #re-add sign
}


# ------------------------------- numformat -----------------------------------#
# source: https://stackoverflow.com/questions/12643391/how-to-remove-leading-0-in-a-numeric-r-variable
# Author: Stefan

numformat <- function(x, digits = 2) { 
  ncode <- paste0("%.", digits, "f")
  sub("^(-?)0.", "\\1.", sprintf(ncode, x))
}

# it displays it in quotes (maybe add noquote())

# ------------------------------- modified TIRT_mplus_syntax -----------------------------------#
# Author: Susanne Frick
# some modifications added by Rebekka Kupffer

#' tirt.mplus.syntax.mod
#'
#' Function to create Mplus syntax for TIRT (only for blocks of 3 items, full ranking,
#' assumes items are ordered blockwise)
#'
#' @param design.load loading matrix of items, rows=items, cols=traits
#' @param names.pairs if NULL they will be named
#' @param item.short shortening to indicate from which questionnaire the items belong
#' @param id.var ID-variable or case-number in quotes
#' @param file.inp file name (directory and name) for Mplus syntax (saved as R-object)
#' @param file.data file name for the data file (.R) without directory: assumed in the same folder as file.inp
#' @param title title for the model
#' @param out.command character string passed to the output command (for example "sampstat standardized;")
#' @param fscores.file name of the file where (and with what name) to save the factor-scores  (e.g., "fscores.dat")
#' @param missings.as how missings are coded in the dataset (e.g., "-99")
#'
#'
#' @return Mplus syntax as list of character vectors, commands as list entries
#' syntax written to file.inp
#' syntax printed to console (with cat)
#' create syntax: syntax <- tirt.mplus.syntax(...)
#' subsequently: cat(paste(syntax, collapse="\n\n"), file="mplus.inp")
#'
#' @export
#'
tirt.mplus.syntax.mod <- function(design.load, names.pairs, item.short, id.var, file.inp=NULL, 
                                  file.data, title, out.command, fscores.file, missings.as,
                                  start.vals.traits=NULL, # starting values for correlations between traits
                                  all.var=NULL, #Variablennamen im Datensatz
                                  use.obs=NULL, #sollen alle observationen verwendet werden, oder nur ein Teil der SP?
                                  add.var=NULL, # Variable, die man zusätzlich noch in die Analysen aufnehmen möchte (muss im Datensatz nach ID und vor den binären Vergleichen kommen)
                                  continuous.var=FALSE, # ist die hinzugefügten Variable continuous?
                                  add.analyses=NULL){ # Analysen, die man noch zusätzlich berechnen möchte
  
  #syntax elements are first created seperately, converted to characters, name: ch.(...)
  #lines are seperated by \n (= new line)
  #use cat() to view result with line breaks
  
  if(is.null(rownames(design.load))) {
    #create itemnames i01, i02 ...
    names.items <- ifelse((1:nrow(design.load))<10,paste0(item.short, "0",1:nrow(design.load)),paste0(item.short,1:nrow(design.load)))
  } else {
    names.items <- rownames(design.load)
  }
  
  
  #design.mat: design matrix of mfc: rows=pairwise comparisons (p), cols=items (i)
  # items assumed to be ordered blockwise:
  # design.block:
  #  i1  i2  i3
  #p1 1   -1  0
  #p2 1   0   -1
  #p3 0   1   -1
  # repeated for each block
  
  design.mat <- matrix(rep(0,nrow(design.load)^2),
                       nrow=nrow(design.load),ncol=nrow(design.load))
  nblocks <- nrow(design.load)/3
  #design.block: part of design matrix for one block
  design.block <- matrix(c(1,-1,0,1,0,-1,0,1,-1),
                         nrow=3,ncol=3,byrow=TRUE)
  #fill blocks in design.mat with block design
  for (i in 1:nblocks) {
    design.mat[(3*i-2):(3*i),(3*i-2):(3*i)] <- design.block
  }
  
  #create pair names if argument empty
  if (is.null(names.pairs)) {
    #create names.pairs: i01i02, i01i03, i02i03 ...
    #for each pairwise comparison: extract itemnames involved, paste together
    #store as rownames(design.mat)
    names.pairs <- NULL
    for (i in 1:nrow(design.mat)) {
      p.name <- names.items[design.mat[i,]!=0]
      names.pairs[i] <- paste0(p.name, collapse="")
    }
  }
  
  #create trait names if colnames(design.load) empty
  if (is.null(colnames(design.load))) {
    names.traits <- paste0("Trait",1:ncol(design.load))
  } else {
    names.traits <- colnames(design.load)
  }
  names.traits
  
  rownames(design.mat) <- names.pairs
  colnames(design.mat) <- names.items
  
  #for variable command: pair names pasted to character
  ch.names.pairs <- paste(names.pairs, collapse="\n")
  ch.all.var <- paste(all.var, collapse="\n")
  ch.add.var <- paste(add.var, collapse="\n")
  
  ######### create model #############
  
  #comp.mat: loadings: rows=comparisons, cols=traits
  comp.mat <- design.mat %*% design.load
  comp.mat.abs <- design.mat %*% abs(design.load)
  comp.mat.abs
  
  ##loadings for pairwise comparisons
  #trait.loadings: character with entries for each trait
  trait.loadings <- NULL
  for (f in 1:ncol(comp.mat)) {
    #comp.trait: part of comparison matrix for trait f
    comp.trait <- comp.mat[grep(1,abs(comp.mat[,f])),f]
    #comp.trait.abs: loading 1/-1 for 1st/2nd item in comparison
    comp.trait.abs <- comp.mat.abs[grep(1,abs(comp.mat.abs[,f])),f]
    
    for (p in 1:length(comp.trait)) {
      #extract pair name
      pair <- names(comp.trait)[p]
      #extract pair loading
      pair.load <- comp.trait[p]
      #extract item name from design.mat
      #select line: rowname=pair
      line.design.mat <- design.mat[grep(pair,rownames(design.mat)),]
      #select item
      item <- names(line.design.mat[line.design.mat==comp.trait.abs[p]])
      #create loading name, e.g. L01, resp. L01_n for items 2nd in comparison
      item.load <- ifelse(comp.trait.abs[p]==1, paste0("L_",item),
                          paste0("L_",item,"_n"))
      
      if (p==1) {
        #first pair loading on trait: begin with 'Trait name BY'
        trait.loadings[f] <- paste0(names.traits[f]," BY \n",pair,"*",pair.load," (",item.load,") \n")
      } else if (p==length(comp.trait)) {
        #last pair loading on trait: end ;
        trait.loadings[f] <- paste0(trait.loadings[f],paste0(pair,"*",pair.load," (",item.load,"); \n"))
      } else {
        trait.loadings[f] <- paste0(trait.loadings[f],paste0(pair,"*",pair.load," (",item.load,") \n"))
      }
    }
  }
  trait.loadings
  
  #as 1 character string (pasted for traits)
  ch.trait.loadings <- paste0(trait.loadings, collapse="\n")
  
  #means for all traits are set to 0
  ch.m.traits <- paste0(paste0("[",names.traits,"@0];", collapse=" "),"\n")
  
  #variances for all traits are set to 1
  ch.var.traits <- paste0(paste0(names.traits,"@1;", collapse=" "),"\n")
  
  #starting values for correlations between traits
  if(is.null(start.vals.traits)){
    trait.cor <- NULL
    if (ncol(design.load)>1) {
      trait.pairs <- combn(names.traits, 2)
      for (i in 1:ncol(trait.pairs)) {
        trait.cor[i] <- paste0(trait.pairs[1,i]," WITH ",trait.pairs[2,i],"*0", collapse=" ")
      }
      trait.cor
      ch.trait.cor <- paste0(paste0(trait.cor,collapse=";\n"),";\n")
    }
  }
  else{
    trait.cor <- NULL
    if (ncol(design.load)>1) {
      trait.pairs <- combn(names.traits, 2)
      for (i in 1:ncol(trait.pairs)) {
        trait.cor[i] <- paste0(trait.pairs[1,i]," WITH ",trait.pairs[2,i],"*", start.vals.traits[i], collapse=" ")
      }
      trait.cor
      ch.trait.cor <- paste0(paste0(trait.cor,collapse=";\n"),";\n")
    }
  }
  
  ##uniquenesses
  #declare uniquenesses and set their starting values (at 2)
  #names uni, e.g. e0102
  names.uni <- gsub(item.short, paste0("e",item.short), names.pairs)
  #names.uni <- sub(item.short,"e",names.uni)
  ch.uni.start <- paste0(names.pairs,"*2 (", names.uni, "); \n", collapse="")
  
  #uniqueness names
  names.uni.items <- gsub(item.short, paste0("e",item.short), names.items)
  
  #declare correlated uniquenesses and set their starting values
  #cor.uni.start: matrix with elements seperately
  cor.uni.start <- matrix(nrow=nrow(design.load),ncol=4)
  paste.cor.uni.start <- NULL
  comb.pairs <- combn(1:3,2)
  for (p in 1:(nrow(design.load)/3)) {
    for (j in 1:3) {
      cor.uni.start[3*p-3+j,1] <- names.pairs[3*p-3+comb.pairs[1,j]] #1st pair
      cor.uni.start[3*p-3+j,2] <- names.pairs[3*p-3+comb.pairs[2,j]] #2nd pair
      cor.uni.start[3*p-3+j,3] <- ifelse(j==2,-1,1) #negative starting value for 2nd pair, see 2011, p.490, (13)
      cor.uni.start[3*p-3+j,4] <- ifelse(j==2,paste0("(",names.uni.items[3*p-3+j],"_n)"),paste0("(",names.uni.items[3*p-3+j],")")) # _n for 2nd pair
      #columns of cor.uni.start pasted together for each line
      paste.cor.uni.start[3*p-3+j] <- paste0(cor.uni.start[3*p-3+j,1]," WITH ",cor.uni.start[3*p-3+j,2],"*",
                                             cor.uni.start[3*p-3+j,3]," ",cor.uni.start[3*p-3+j,4],"; \n")
    }
  }
  cor.uni.start
  #paste all lines together
  ch.cor.uni.start <- paste0(paste.cor.uni.start,collapse="")
  
  ##contraints
  #factor loadings relating to the same item are equal in absolute value
  #for every 2nd item in a triplet
  second.items <- names.items[seq(2,nrow(design.load),by=3)]
  neg.loads <- paste(paste0(sub(item.short, paste0("L_",item.short),second.items),"_n = ",
                            paste0("-",sub(item.short, paste0("L_",item.short), second.items)), "; \n"),collapse="")
  
  ch.neg.loads <- paste0("MODEL CONSTRAINT: \n\n! factor loadings relating to the same item are equal in absolute magnitude \n",neg.loads)
  
  
  #pairs uniqueness is equal to sum ot 2 utility uniquenesses
  #matrix with elements seperately
  uni.pair.util <- matrix(nrow=nrow(design.mat),ncol=3)
  # pair uniquenesses
  uni.pair.util[,1] <- paste(names.uni,"= ",sep=" ")
  # utility uniquenesses (without parantheses)
  uni.utils <- sub("[(]","",cor.uni.start[,4]) # [] around strings with logical meaning in R
  uni.utils <- sub(")","",uni.utils)
  
  for (b in 1:nblocks) {
    for (j in 1:3) {
      # 2 utility uniquenesses involved in pair, - for 2nd
      uni.pair.util[3*b-3+j,2] <- ifelse(comb.pairs[1,j]==2,
                                         paste0("-",uni.utils[3*b-3+comb.pairs[1,j]]),
                                         uni.utils[3*b-3+comb.pairs[1,j]])
      uni.pair.util[3*b-3+j,3] <- ifelse(comb.pairs[2,j]==2,
                                         paste0(" - ",uni.utils[3*b-3+comb.pairs[2,j]],"; \n"),
                                         paste0(" + ", uni.utils[3*b-3+comb.pairs[2,j]],"; \n"))
    }
  }
  ch.uni.pair.util <- paste0(t(uni.pair.util), collapse = "")
  
  #fix one uniquness per block for identification
  ch.uni.util.fix <- paste0(uni.utils[seq(1,nrow(design.mat),3)],"=1; \n",collapse="")
  
  #------------------------------------------------------------------------------------#
  
  #combine to model command
  if (ncol(design.load)>1) {
    tirt.model <- paste(ch.trait.loadings, "! means for all traits are set to 0",ch.m.traits,
                        "! variances for all traits are set to 1",ch.var.traits,
                        "! starting values for correlations between traits",ch.trait.cor,
                        "! declare uniquenesses and set their starting values",ch.uni.start,
                        "! declare correlated uniquenesses and set their starting values",ch.cor.uni.start,
                        sep="\n")
  } else {
    tirt.model <- paste(ch.trait.loadings, "! means for all traits are set to 0",ch.m.traits,
                        "! variances for all traits are set to 1",ch.var.traits,
                        "! declare uniquenesses and set their starting values",ch.uni.start,
                        "! declare correlated uniquenesses and set their starting values",ch.cor.uni.start,
                        sep="\n")
  }
  tirt.model.constraint <- paste(ch.neg.loads,
                                 "! pair's uniqueness is equal to the sum of 2 utility uniquenesses",ch.uni.pair.util,
                                 "! fix one uniquness per block for identification",ch.uni.util.fix,
                                 sep="\n")
  
  ##################################
  categorial <- ifelse(continuous.var == TRUE, 
                       paste0("\nCATEGORICAL ARE ", names.pairs[1],"-",names.pairs[length(names.pairs)], ";"),
                       "\nCATEGORICAL ARE ALL;")
  
  observations <- if(!is.null(use.obs)) paste0("\nUSEOBSERVATIONS = ", use.obs, ";")
  
  ##################################
  
  #combine syntax elements to list
  syntax.tirt <-  list(paste0("TITLE: ",title,";"),
                       paste0("DATA: FILE IS '",file.data,"';"),
                       paste0("VARIABLE: \nNames ARE \n", id.var, "\n", ch.all.var, "\n", ch.names.pairs,";\n",
                              "USEVARIABLES ARE \n", ch.add.var, "\n", ch.names.pairs,";",
                              categorial, observations, "\nIDVARIABLE = ", id.var, ";\n",
                              "MISSING ARE ALL (", missings.as, ");\n"),
                       paste0("ANALYSIS: \nESTIMATOR = ULSMV;\nPARAMETERIZATION = theta;"),
                       paste0("MODEL: \n",tirt.model, "\n", add.analyses),
                       tirt.model.constraint,
                       paste0("OUTPUT:\n",out.command),
                       paste0("SAVEDATA:\n", "FILE IS '", fscores.file, "';\n", "SAVE = FSCORES;"))
  
  names(syntax.tirt) <- c("title","data","variable","analysis","model","model.constraint","output", "savedata")
  
  #save syntax list as R object
  #can be loaded in later session to manipulate list entries
  if(is.null(file.inp)==FALSE) save(syntax.tirt, file=file.inp)
  #print syntax to console
  cat(paste(syntax.tirt, collapse="\n\n"))
  #return syntax list
  return(syntax.tirt)
}

