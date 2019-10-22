dissimilarity <- function(col1,col2){
  
  #Total pops
  col1sum = sum(col1)
  col2sum = sum(col2)
  
  #https://en.wikipedia.org/wiki/Index_of_dissimilarity
  0.5 * (sapply(1:length(col1), function(x) {
    abs((col1[x]/col1sum) - (col2[x]/col2sum))
  })) %>% unlist %>% sum
  
}
