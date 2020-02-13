library(spdep)

#Adapted from Nema's ACD code in permutation tests
aacd <- function(data.sf,var,expz = 1){
  
  names(data.sf)[names(data.sf)==var] = 'var'  
  
  #Get absolute contig diff values (via Nema code in permutation_tests)
  data.sp <- as_Spatial(data.sf)
  neighbours <- poly2nb(data.sp, queen=F)
  
  data.n.mx <- nb2mat(neighbours,style = 'B')
  
  x.diff.matrix <- outer(data.sp$var,data.sp$var,FUN="-")
  x.diff.matrix <- abs(x.diff.matrix)
  x.diff.matrix <- x.diff.matrix^expz#Raise result to power
  
  #To avoid having to find unique values (which won't work here, too many similar)
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



#values: vector of values for a geographical space
#neighbours: neighbour list for that space via spdep
weightedAACD = function(values, neighbours){


#CRASH
#This was the code I wrote for tracking down the crash issue...
#Do same for loops in R, look for issues
#This is the key line:
#all we're looking for: j neighbours of neighbours (2nd order) that are also i neighbours (first order).
#(NOT overlapping neighbours of neighbours of j and neighbours of neighbours of i - they don't overlap!)

#This is the version that just demonstrates we can find the right grid squares...
#Count border calcs we carry out
bordercount <<- 0
#AACD or some other measure. Sum, then find mean for all borders
bordercalc <<- 0

#temp: check actual values, see if order matches C++. Should, in theory...
valz <- list()

for(i in 1:length(values)){
  
  
  #Get neighbours of cell i
  i_nb <- neighbours[[i]]
  
  
  #For each neigbour of cell i: j
  for(j in 1:length(i_nb)){
    
    
    #Avoid double-running e.g. i, j and j, i
    if(i < i_nb[[j]]){
      
      
      bordercount = bordercount + 1
      
      #Get neighbours of cell j
      j_nb <- neighbours[[ i_nb[[j]] ]]
      
      #We then want neighbours of neighbours of cell j
      #all we're looking for: j neighbours of neighbours (2nd order) that are also i neighbours (first order).
      #We're sticking to rook contiguity - that'll make more sense for real geogs
      #- so those diagonals that are neighbours of i are not immediate neighs of j
      j_nb_of_nb <- lapply(j_nb, function(x) neighbours[[x]])
      names(j_nb_of_nb) <- j_nb
      
      #Need to exclude i's neighbours, otherwise neigh of neigh of j gets opposite-side i neighbour we don't want
      #Neighbours of neighbours of j includes i's neighbours.
      #We'll pick up the ones we need from from the other two neighbours of neighbours of j.
      j_nb_of_nb[[paste0(i)]] <- NULL
      
      x <- unlist(j_nb_of_nb)
      
      #remove extra digit that unlist adds so we get j neighbour name properly
      names(x) <- str_sub(names(x), end=-2)
      
      
      #What's in both?
      #If we exclude j itself, I think this works
      ij_bordershared <- i_nb[which(i_nb %in% x)]
      ij_bordershared <- ij_bordershared[which(ij_bordershared!=i_nb[[j]])]
      
      #And then get the attached names for those two...
      pairs <- as.numeric(names(x)[which(x %in% ij_bordershared)])
      
      
      
      #So: i and j are our main pair
      #Then we have kl and mn. Find ACD for each of those.
      ij_acd <- abs(values[i] - values[i_nb[[j]]])
      
      # cat('i, j: ', i, ',',i_nb[[j]],'\n')
      # cat('k, l: ', ij_bordershared[1], ',',pairs[1],'\n')
      
      
      klmn_acd = list()
      klmn_acd[[1]] <- abs(values[ ij_bordershared[1] ] - values[  pairs[1] ])
      
      #We may only have one extra pair if i is on edge...
      if(length(ij_bordershared)==2){
        klmn_acd[[2]] <- abs(values[ ij_bordershared[2] ] - values[  pairs[2] ])
        # cat('m, n: ', ij_bordershared[2], ',',pairs[2],'\n')
        
      }
      
      # cat('\n')
      
      
      # cat('ij acd: ', ij_acd, '\n')
      # 
      # for(i in klmn_acd){
      # cat('klmn acd: ',  i, '\n')
      # }
      # 
      
      #Original Nema weighted border calculation:
      #ij are pair. mn and kl are the neighbouring pairs on each side
      #If both mn kl ACD are smaller than ij ACD
      #Then swap for larger of mn kl
      #otherwise keep as is
      if(max(unlist(klmn_acd)) < ij_acd){
        
        final_acd = max(unlist(klmn_acd))
        
      } else {
        
        final_acd = ij_acd
        
      }
      
      bordercalc = bordercalc + final_acd
      
      valz[[length(valz)+1]] <- final_acd
      
    }
    
  }
  
}

  #final mean border calc, av for all borders (halved)
  bordercalc = bordercalc / bordercount
  
  return(list(bordercalc,unlist(valz)))

}






#Local Morans I on each side of border
#vector of values from the geography
#(Rook) neighbour list to find bordering cells/zones
#Queen neighbour list to put together little surrounding areas more easily
#Original spatial dataset (in sp-friendly format) to get those surrounding areas
localMoransI_onEachSideOfBorder = function(values,neighbours,neighbours.queen,origSpatial){
  
  bordercount <<- 0
  #AACD or some other measure. Sum, then find mean for all borders 
  bordercalc <<- 0 
  
  #temp: check actual values, see if order matches C++. Should, in theory...
  valz <- list()
  
  for(i in 1:length(values)){
    
    # cat("i:", i, "\n")
    
    
    #Get neighbours of cell i
    i_nb <- neighbours[[i]]
    
    #Get queen-contig neighbours of cell i  
    i_nb.q <- neighbours.queen[[i]]
    
    
    #For each neigbour of cell i: j
    for(j in 1:length(i_nb)){ 
       
      # cat("j:", j, "\n")
      
      #Avoid double-running e.g. i, j and j, i
      if(i < i_nb[[j]]){
        
        bordercount = bordercount + 1
        
        #We're finding two Moran's Is
        #For the nine cells of i and 8 surrounding 
        #And nine cells of j and 8 surrounding 
        
        i_local <- origSpatial[c(i,i_nb.q),] 
        j_local <- origSpatial[c(i_nb[[j]], neighbours.queen[[  i_nb[[j]]  ]] ),]
        
        i_local.nbs <- poly2nb(i_local,queen=F)
        j_local.nbs <- poly2nb(j_local,queen=F)
        
        #This is relying on these numbers being the same in the spatial thingyo... watch out
        #Easiest way though.
        
        #Check if local cell values are all identical. If they are, no sensible Moran value (function throws error)
        #Just set to zero - doesn't affect overall score either way in that case.
        # cat("i:", i,", j:",j,"\n")
        
        i_moran = -999
        j_moran = -999
        
        
        if(sd(i_local$peeps1)>0){
          i_moran <- as.numeric(moran.test(i_local$peeps1, nb2listw(i_local.nbs))$estimate[1])
        } else{
          i_moran <- 0
        }
        
        if(sd(j_local$peeps1)>0){
          j_moran <- as.numeric(moran.test(j_local$peeps1, nb2listw(j_local.nbs))$estimate[1])
        } else{
          j_moran <- 0
        }
          


        #Sum those two for final border score
        #Also test make negative so that e.g. smooth Moran results (that would be positive) push in same direction as
        #Low border differences / ACDs
        # bordercalc <- bordercalc + -(i_moran + j_moran)
        bordercalc <- bordercalc + i_moran + j_moran
        
        # valz[[length(valz)+1]] <- -(i_moran + j_moran)
        valz[[length(valz)+1]] <- i_moran + j_moran
        
      }#end if
      
    }#end for j
    
  }#end for i
  
  return(list(mean = bordercalc/bordercount,values = unlist(valz)))
  
}#end of the world





#sf spatial data containing value and grid geography
#varname to use (is here assuming just one var and var2 is 1-var...)
centralisationIndex_onEachSideOfBorder = function(data.sf,var){
  
  names(data.sf)[names(data.sf)==var] <- "var"
  
  data.sp <- as_Spatial(data.sf)
  
  #Rook contig plz! We only want bordering cells
  neighbours <- poly2nb(data.sp, queen=F)
  #Except we also want to get little group surrounding cell...
  neighbours.queen <- poly2nb(data.sp, queen=T)
  
  
  bordercount <<- 0
  #AACD or some other measure. Sum, then find mean for all borders
  bordercalc <<- 0 
  
  #Get each border/neighbour pair value
  valz <- list()
  
  
  for(i in 1:nrow(data.sf)){
    
    # cat("i:", i, "\n")
    
    
    #Get neighbours of cell i
    i_nb <- neighbours[[i]]
    
    #Get queen-contig neighbours of cell i 
    i_nb.q <- neighbours.queen[[i]]
    
    
    #For each neigbour of cell i: j
    for(j in 1:length(i_nb)){  
       
      # cat("j:", j, "\n")
      
      #Avoid double-running e.g. i, j and j, i
      if(i < i_nb[[j]]){
        
        bordercount = bordercount + 1
        
        #We're finding two lots of RCI
        #For the nine cells of i and 8 surrounding
        #And nine cells of j and 8 surrounding
        
        i_local <- data.sf[c(i,i_nb.q),] 
        j_local <- data.sf[c(i_nb[[j]], neighbours.queen[[  i_nb[[j]]  ]] ),]
        
        #need to use centroids for distances
        i_local.centroids <- st_centroid(i_local)
        j_local.centroids <- st_centroid(j_local)
        
        #Distance from central point, which is indexed either by i or j
        #Can use ID field to identify
        i_local$distances <- st_distance(i_local.centroids[i_local.centroids$id==i,],i_local.centroids,by_element = T)
        j_local$distances <- st_distance(j_local.centroids[j_local.centroids$id==i_nb[[j]],],j_local.centroids,by_element = T)
        
        #They both look like sensible distances
        #Cos gridsquares, will be some the same (three the same if group of 9)
        #Centre, rook and queen cells.
        
        #Should be enough to get rci with Meng Le's function
        #Recall: y should be central for it to be positive
        #So gonna reverse to keep sensible-ish
        j_rci <- rci(x = 1-j_local$var,y=j_local$var,sort.var = j_local$distances)
        i_rci <- rci(x = 1-i_local$var,y=i_local$var,sort.var = i_local$distances)
        
        # cat('i_rci: ',i_rci,'\n')
        # cat('j_rci: ',j_rci,'\n')
        # 
        #Final result: absolute value of subtracting both
        
        #Sum those two for final border score
        bordercalc <- bordercalc + abs(i_rci-j_rci)
        
        valz[[length(valz)+1]] <- abs(i_rci-j_rci)
        
      }#end if
      
    }#end for j
    
  }#end for i
  
  return(list(mean=bordercalc/bordercount,values=unlist(valz)))
  
}#end of the world





#ASSUMES ALREADY ORDERED BY DISTANCE
relativeCentralisationIndex = function(data, varxname, varyname){
  
  #Straight strings, swap
  names(data)[names(data)==varxname] = 'varx'
  names(data)[names(data)==varyname] = 'vary'
  
  tot = list()
  
  #Sum cum proportions multiplied for first lot
  for(i in 2:nrow(data)){
    
    x = cumsum(data[1:(i-1),'varx'])
    x = x[length(x)]#get final cumulative sum
    
    y = cumsum(data[1:(i),'vary'])
    y = y[length(y)]#get final cumulative sum
    
    tot[[length(tot)+1]] <- x*y
    
  }
  
  #Sum all those results
  tot = sum(unlist(tot))
  
  tot2 = list()
  
  #Sum cum proportions multiplied for first lot
  for(i in 2:nrow(data)){
    
    x = cumsum(data[1:(i),'varx'])
    x = x[length(x)]#get final cumulative sum
    
    y = cumsum(data[1:(i-1),'vary'])
    y = y[length(y)]#get final cumulative sum
    
    tot2[[length(tot2)+1]] <- x*y
    
  }
  
  #Sum all those results
  tot2 = sum(unlist(tot2))
  
  tot-tot2
  
}




#MENG LE'S 25* FASTER RCI FUNCTION
#BE AWARE ORDER IS REVERSED I.E. IF Y IS CENTRALISED, RCI WILL BE POSITIVE
rci <- function(x, y, sort.var) {
  ordering <- order(sort.var)
  N <- length(ordering)
  x <- x[ordering]
  y <- y[ordering]
  
  ##  cumulative sums
  csum.x <- cumsum(x) / sum(x)
  csum.y <- cumsum(y) / sum(y)
  
  ## The below has been checked to be right;
  out <- t(csum.y[-N]) %*% csum.x[-1] - t(csum.y[-1]) %*% csum.x[-N]
  return(out)
}







avInVicinityForEachBorderPair <- function(queen.matrix, var){
  
  #Check it's row-normalised. Tick.
  #apply(m,1,sum)
  spatial.lag <- queen.matrix %*% var
  
  #Spatial lag is for each cell in order
  #But need to have that value for each ij pair, just returning for i
  #Easiest way probably to run through usual order of neighbours
  lag4calc = list()
  for(i in 1:length(var)){
    for(j in 1:length(neighbours[[i]])){
      if(i< neighbours[[i]][j]){  
        lag4calc[[length(lag4calc)+1]] <- spatial.lag[i]
      }
    }
  }
  
  unlist(lag4calc)
  
}




acd_variogram <- function(grid,varname){
  
  names(grid)[names(grid)==varname] = "var"
  
  dists <- st_distance(st_centroid(grid),st_centroid(grid), by_element = F)
  
  #I'm gonna assume it's doing x1, all y, x2, all y...
  #That is looking correct.
  dists.data <- data.frame(
    x=grid$var[rep(1:nrow(grid),each=nrow(grid))],
    y=grid$var[rep(1:nrow(grid),times=nrow(grid))],
    distance = as.numeric(dists)
  )
  
  #Absolute diff
  dists.data <- dists.data %>% 
    mutate(diff = abs(x-y))
  
  #Now we can bundle all our 65536 abs diffs into distance bands
  #Actually, what's the distr of the dists?
  # plot(density(dists.data$distance))
  # length(unique(dists.data$distance))
  # 
  # #And the distr of the diffs? (Not at all confusing language there!)
  # plot(density(dists.data$diff))
  # hist(dists.data$diff)
  # 
  #Add distance bands for examining groups at different distances
  #Min distance is 500 so this works as band is 530. Could use cut_width directly too.
  dists.data <- dists.data %>% 
    mutate(distgroup = cut_interval(distance,20))
  
  #Can now group and do whatever. So let's try a few things
  dists.grouped <- dists.data %>% 
    group_by(distgroup) %>% 
    summarise(mean = mean(diff), max = max(diff), min = min(diff), sd = sd(diff))
  
  ggplot(dists.grouped %>% gather(key=summaryvar,value=value,mean:sd),
         aes(x = distgroup, y = value, colour= summaryvar, group=summaryvar)) +
    geom_point() +
    geom_line()
  
}










