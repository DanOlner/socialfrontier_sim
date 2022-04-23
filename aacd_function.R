aacd <- function(data.sf,var,expz = 1){
  
  names(data.sf)[names(data.sf)==var] = 'var'  
  
  #Get absolute contig diff values (via Nema code in permutation_tests)
  data.sp <- as_Spatial(data.sf)
  neighbours <- poly2nb(data.sp, queen=F)
  
  data.n.mx <- nb2mat(neighbours,style = 'B')
  
  x.diff.matrix <- outer(data.sp$var,data.sp$var,FUN="-")
  x.diff.matrix <- abs(x.diff.matrix)
  x.diff.matrix <- x.diff.matrix^expz#Raise result to power
  
  #To avoid having to find unique values (which won't work here I don't think, too many similar)
  #Keep only upper of neighbours
  data.n.mx[lower.tri(data.n.mx)] <- 0
  
  ## paired differences amongst geographical neighbours
  #Which is the right number - 97 borders. Good good.
  x.diff.W0 <- x.diff.matrix[data.n.mx==1]
  #x.diff.W0 <- unique(x.diff.W0)
  
  #AACD is then just that value divided by border number. Or just... mean.
  mean(x.diff.W0)
  
  return(list(mean = mean(x.diff.W0), values = x.diff.W0))
  
}