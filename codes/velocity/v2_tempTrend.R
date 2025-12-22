# Slope Function for Linear Model for Time-Series Data
slpFUN <- function(k) {  
  timestamps <- 1:nlayers(k) 
  mn <- calc(k,
             fun = function(x) {
               if (all(is.na(x)))
                 return(NA)
               else
                 return(coef(lm(x ~ timestamps))[2])})
  return(mn <- (mn*length(timestamps))/30)
}
