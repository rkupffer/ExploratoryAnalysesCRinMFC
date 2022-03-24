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

