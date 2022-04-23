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
#By default, just first order lag for the local RCI measure thingyo.
centralisationIndex_onEachSideOfBorder = function(data.sf,var,lag=1){
  
  names(data.sf)[names(data.sf)==var] <- "var"
  
  data.sp <- as_Spatial(data.sf)
  
  #Rook contig plz! We only want bordering cells
  neighbours <- poly2nb(data.sp, queen=F)
  
  #This is an index for the local areas we want to use in the centralisation index calcs
  neighbours.queen <- poly2nb(data.sp, queen=T)
  
  if(lag>1){
    
    x <- nblag(neighbours.queen,maxlag = lag)
    neighbours.queen <- nblag_cumul(x)#Combine into one list
    
  }
  
  # cat('typical neighbours.queen size: ',length(neighbours.queen[[150]]),'\n')
  
  #Tick
  # data.sf$check=0
  # data.sf$check[neighbours.queen[[150]]]=1
  # plot(data.sf)
  
  bordercount <<- 0
  #AACD or some other measure. Sum, then find mean for all borders
  bordercalc <<- 0 
  
  #Get each border/neighbour pair value
  valz <- list()
  
  
  for(i in 1:nrow(data.sf)){
    
    # cat("i:", i, "\n")
    
    #Get neighbours of cell i
    i_nb <- neighbours[[i]]
    
    #Get queen-contig neighbours (or higher order list of neighbours of neighbours of...) of cell i 
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






#This is doing what we want at the moment - 
#Finding proportion var in area -
#Only cos proportions are used as the var. Would break with actual counts.
avInVicinityForEachBorderPair <- function(queen.matrix, neighbours, var, exponent=1){
  
  #Check it's row-normalised. Tick.
  #apply(queen.matrix,1,sum)
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
  
  lag4calc <- unlist(lag4calc)
  lag4calc <- lag4calc^exponent
  
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
    geom_line() +
    theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
}







#Adjacency matrix neighbour list, including option to only return unique pairs
neighbourlist <- function(sf_object, queen = F, uniquepairs = F){
  
  if(uniquepairs){
    obj.spdep <- as_Spatial(sf_object)
    neighbours <- poly2nb(obj.spdep,queen = queen)
    matrx <- nb2mat(neighbours,style = 'B')
    matrx[lower.tri(matrx)] <- 0
    x <- mat2listw(matrx)
    return(x[[2]])
    
  } else {
    
    obj.spdep <- as_Spatial(sf_object)
    return(poly2nb(obj.spdep,queen = queen))
    
  }
  
}


#WAS USING THIS TO ITERATE OVER CUTOFF VALUES - USE GENERIC FRONTIER LENGTH FUNCTION BELOW INSTEAD
#Connected frontier networks for e.g. finding frontier lengths
#Currently hardcoded to the peeps1 variable in the incoming sf dataframe
#Returns an igraph
# frontiernetworks <- function(frontiercutoff, spatial.sf){
# 
# 
#   neighbs.spdep <- as_Spatial(spatial.sf)
#   neigbs.rook <- poly2nb(neighbs.spdep,queen = F)
#   neigbs.queen <- poly2nb(neighbs.spdep,queen = T)
#   
#   #Let's test getting all relevant connected frontiers, assuming (5,6) has a frontier
#   #What queen contig neighbours do 5 and 6 share?
#   #1,2,9,10. Yeah, they're all the potential linked frontier zones.
#   
#   #Not making generic as regards the data column for frontier-finding yet
#   connectedfrontiers_edgelist = list()
#   #Make sure to empty this before starting
#   # potentiallink_pairlist <- list()
#   
#   #Loop 1: across each main zone...
#   for(mainzone in 1:nrow(gridedge_centre)){
#     # for(mainzone in 88){
#     
#     #And each zone pair made from that zone's rook neighbours
#     for(mainzone_neighbour in neigbs.rook[[mainzone]]){
#       # for(mainzone_neighbour in 104){
#       
#       # cat('mainzone: ',mainzone, ", mainzone_neighbour: ", mainzone_neighbour, "\n")
#       
#       #Is that top-level pair a frontier (based on our def)?
#       #If not, can skip this check for connected frontiers
#       if(abs(gridedge_centre$peeps1[mainzone] - gridedge_centre$peeps1[mainzone_neighbour]) > frontiercutoff){
#         
#         # cat('Past if frontier test -- mainzone: ',mainzone, ", mainzone_neighbour: ", mainzone_neighbour, "\n")
#         
#         
#         #Then check around that pair for any of the six possible connecting frontiers
#         #(Or possibly more if it's irregular polygons; I *think* the same code will work for both but need to test)
#         #If there are connecting frontiers, add to connected frontier edge list
#         
#         #Get all the potential linked frontier zones.
#         bothqueen <- neigbs.queen[[mainzone]][neigbs.queen[[mainzone]] %in% neigbs.queen[[mainzone_neighbour]]]
#         
#         #Include two main pair zones
#         bothqueen <- c(bothqueen,mainzone,mainzone_neighbour)
#         
#         #Cycle through each, extract pairs that have rook contig WITH ANY OTHER IN THE SAME LIST
#         potentiallink_pairlist <- list()
#         
#         for(i in bothqueen){
#           
#           rookoverlap <- bothqueen[bothqueen %in% neigbs.rook[[i]]]
#           #Exclude the original ones we're checking
#           rookoverlap <- rookoverlap[!rookoverlap %in% c(mainzone,mainzone_neighbour)]
#           
#           #Now we can add any we found to the potentiallink_pairlist
#           for(j in rookoverlap){
#             
#             #check if there's a frontier and add if so
#             if(abs(gridedge_centre$peeps1[j] - gridedge_centre$peeps1[i]) > frontiercutoff){
#               
#               #Keep same order so duplicates can be removed
#               if(j > i){
#                 potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(i,j)
#               } else {
#                 potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(j,i)
#               }
#               
#               
#             }#end if frontiercutoff
#             
#           }#end for j
#           
#         }#end for i
#         
#         #Keep only unique pairs
#         potentiallink_pairlist <- unique(potentiallink_pairlist)
#         
#       } else {#Else if this main / neighbour pair are NOT frontier...
#         
#         #Hack to empty list if not a frontier
#         #otherwise if a frontier on prev. iteration, it keeps prev. list items
#         #If the logic here were better, this wouldn't be necessary!
#         potentiallink_pairlist <- list()
#         
#       }
#       
#       # print(paste0("potentiallink_pairlist length:",length(potentiallink_pairlist)))
#       
#       if(length(potentiallink_pairlist) > 0){
#         
#         for(k in potentiallink_pairlist){
#           
#           # cat("k: ",k,"\n")
#           
#           #Keep same order for pairs of pairs, so can keep unique edge links
#           #First, put the original/main pairlist in order. Suspect ifelse won't work... ah, I see why now reading the help
#           #Store each pair as a string
#           if(mainzone < mainzone_neighbour){
#             mainpair <- paste0(mainzone,",",mainzone_neighbour)
#           } else {
#             mainpair <- paste0(mainzone_neighbour,",",mainzone)
#           }
#           
#           #And keep in same order - doesn't matter what method used as long as pairs are in same order
#           #Is not a directed graph, we just need unique pairs
#           if(k[1] + k[2] < mainzone + mainzone_neighbour){
#             
#             # connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <-
#             #   paste0(k[1],",",k[2],"--",mainpair)
#             connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <-
#               c(paste0(k[1],",",k[2]),mainpair)
#             
#           } else {
#             
#             # connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <- 
#             #   paste0(mainpair,"--", k[1],",",k[2])
#             connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <- 
#               c(mainpair,paste0(k[1],",",k[2]))
#             
#           }
#           
#         }
#         
#       }#end if length potentiallink_pairlist
#       
#       
#       
#     }#For mainzone_neighbour (add to edgelist before this)
#     
#   }#For mainzone
#   
#   connectedfrontiers_edgelist <- unique(connectedfrontiers_edgelist)
#   
#   edges.as.matrix <- matrix(unlist(connectedfrontiers_edgelist), ncol = 2, byrow = T)
#   frontier.graph <- igraph::graph_from_edgelist(edges.as.matrix, directed=F)
# 
# }








#Return:
#2 different igraphs - 
#(1) all connected frontiers; 
#(2) All 'consistent' frontiers (frontiers all on same 'side')
#And:
#The frontier shapefile, including a column denoting the component ID and membership count
find.frontier.lengths <- function(frontiercutoff,sourcemap,valuecol){
  
  #One name plz
  names(sourcemap)[names(sourcemap)==valuecol] <- "valuecol"
  
  #1. MAKE FRONTIER SF DATAFRAME FOR VISUALISING THEM
  uniquerookneighbours <- neighbourlist(sourcemap,queen = F,uniquepairs = T)
  frontiers <- list()
  
  #For each zone...
  #Have to do minus here if using unique pairs
  #as the final bottom-right zone has no neighbours after others took those pairs
  for(i in 1:(nrow(sourcemap)-1)){
    
    #And for each neighbour of that zone...
    for(neigh in uniquerookneighbours[[i]]){
      
      #If absolute contig difference is above value for defining as frontier....
      if(abs(sourcemap$valuecol[i] - sourcemap$valuecol[neigh]) > frontiercutoff){
        
        frontiers[[length(frontiers)+1]] <- st_intersection(sourcemap[i,],sourcemap[neigh,])
        
      }
      
    }
    
  }
  
  #https://geocompr.github.io/presentations/attr.html#6 - weird bind_rows error
  frontiers <- do.call(rbind, frontiers)
  
  
  #~~~~~~~~~~~~~~~
  #2. FIND FULL FRONTIER CONNECTION LIST ON 16*16 WITH FRONTIER CUTOFF
  
  #Assemble an edgelist that igraph can turn into a graph with graph_from_edgelist
  #Each edge is a frontier
  #Specified by the two zone numbers on each side e.g. (5,6)
  #Keep same zone number order so each zonepair is a unique edge
  
  #Just use one neighbour list for each of the three same-dimension grids
  neighbs.spdep <- as_Spatial(sourcemap)
  neigbs.rook <- poly2nb(neighbs.spdep,queen = F)
  neigbs.queen <- poly2nb(neighbs.spdep,queen = T)
  
  #Not making generic as regards the data column for frontier-finding yet
  connectedfrontiers_edgelist = list()
  #Make sure to empty this before starting
  
  #Loop 1: across each main zone...
  for(mainzone in 1:nrow(sourcemap)){
    # for(mainzone in 88){
    
    #And each zone pair made from that zone's rook neighbours
    for(mainzone_neighbour in neigbs.rook[[mainzone]]){
      # for(mainzone_neighbour in 104){
      
      cat('mainzone: ',mainzone, ", mainzone_neighbour: ", mainzone_neighbour, "\n")
      
      #Is that top-level pair a frontier (based on our def)?
      #If not, can skip this check for connected frontiers
      if(abs(sourcemap$valuecol[mainzone] - sourcemap$valuecol[mainzone_neighbour]) > frontiercutoff){
        
        cat('Past if frontier test -- mainzone: ',mainzone, ", mainzone_neighbour: ", mainzone_neighbour, "\n")
        
        
        #Then check around that pair for any of the six possible connecting frontiers
        #(Or possibly more if it's irregular polygons; I *think* the same code will work for both but need to test)
        #If there are connecting frontiers, add to connected frontier edge list
        
        #Get all the potential linked frontier zones.
        bothqueen <- neigbs.queen[[mainzone]][neigbs.queen[[mainzone]] %in% neigbs.queen[[mainzone_neighbour]]]
        
        #Include two main pair zones
        bothqueen <- c(bothqueen,mainzone,mainzone_neighbour)
        
        #Cycle through each, extract pairs that have rook contig WITH ANY OTHER IN THE SAME LIST
        potentiallink_pairlist <- list()
        
        for(i in bothqueen){
          
          rookoverlap <- bothqueen[bothqueen %in% neigbs.rook[[i]]]
          #Exclude the original ones we're checking
          rookoverlap <- rookoverlap[!rookoverlap %in% c(mainzone,mainzone_neighbour)]
          
          #Now we can add any we found to the potentiallink_pairlist
          for(j in rookoverlap){
            
            #check if there's a frontier and add if so
            if(abs(sourcemap$valuecol[j] - sourcemap$valuecol[i]) > frontiercutoff){
              
              #Keep same order so duplicates can be removed
              if(j > i){
                potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(i,j)
              } else {
                potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(j,i)
              }
              
              
            }#end if frontiercutoff
            
          }#end for j
          
        }#end for i
        
        #Keep only unique pairs
        potentiallink_pairlist <- unique(potentiallink_pairlist)
        
      } else {#Else if this main / neighbour pair are NOT frontier...
        
        #Hack to empty list if not a frontier
        #otherwise if a frontier on prev. iteration, it keeps prev. list items
        #If the logic here were better, this wouldn't be necessary!
        potentiallink_pairlist <- list()
        
      }
      
      # print(paste0("potentiallink_pairlist length:",length(potentiallink_pairlist)))
      
      if(length(potentiallink_pairlist) > 0){
        
        for(k in potentiallink_pairlist){
          
          # cat("k: ",k,"\n")
          
          #Keep same order for pairs of pairs, so can keep unique edge links
          #First, put the original/main pairlist in order. Suspect ifelse won't work... ah, I see why now reading the help
          #Store each pair as a string
          if(mainzone < mainzone_neighbour){
            mainpair <- paste0(mainzone,",",mainzone_neighbour)
          } else {
            mainpair <- paste0(mainzone_neighbour,",",mainzone)
          }
          
          #And keep in same order - doesn't matter what method used as long as pairs are in same order
          #Is not a directed graph, we just need unique pairs
          if(k[1] + k[2] < mainzone + mainzone_neighbour){
            
            # connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <-
            #   paste0(k[1],",",k[2],"--",mainpair)
            connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <-
              c(paste0(k[1],",",k[2]),mainpair)
            
          } else {
            
            # connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <- 
            #   paste0(mainpair,"--", k[1],",",k[2])
            connectedfrontiers_edgelist[[length(connectedfrontiers_edgelist)+1]] <- 
              c(mainpair,paste0(k[1],",",k[2]))
            
          }
          
        }
        
      }#end if length potentiallink_pairlist
      
      
      
    }#For mainzone_neighbour (add to edgelist before this)
    
  }#For mainzone
  
  connectedfrontiers_edgelist <- unique(connectedfrontiers_edgelist)
  
  #Graph from edge list, find average frontier length
  cat('connectedfrontiersedgelist length ',length(connectedfrontiers_edgelist),'\n')
  
  #Break out of function if we found no frontiers
  if(length(connectedfrontiers_edgelist)==0){
    # return(NULL)
    return(list(frontier.graph = graph.empty(n=0, directed=F), 
                frontier.graph.valid = graph.empty(n=0, directed=F), 
                frontiers.sf = frontiers))
  }
  
  edges.as.matrix <- matrix(unlist(connectedfrontiers_edgelist), ncol = 2, byrow = T)
  frontier.graph <- graph_from_edgelist(edges.as.matrix, directed=F)
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #3. CONSISTENT NEIGHBOURING FRONTIERS----
  
  #That is: using the network already found, and again working in pairs:
  #Check zones on each side of the borders are consistent 
  #So that "this side" of the frontier for each pairs is the same.
  
  #For corner links, needs something different:
  #Use the cutoff, and check that the single corner cell is 
  #On the opposite side of the cutoff from the other 3 cells
  #We already know the rook contig links will be - 
  #but use queen contig to check the corner.
  
  #Those give slightly different info - think I can work with the first
  #https://stackoverflow.com/questions/52159669/how-to-get-the-vertices-of-a-selected-edge-in-igraph-r
  # get.edgelist(frontier.graph)[1,]
  
  #So, we have our two types of pair:
  #Straight line connection
  #Corner connection
  
  #Set entire 'valid edge' attribute; write over with T if valid
  #(Don't think that works... will overwrite below)
  E(frontier.graph)$validedge = F
  
  
  #Cycle through all edges (all pairs of connected frontiers)
  for(i in 1:ecount(frontier.graph)){
    
    E(frontier.graph)[i]$validedege <- F
    
    #We can test for a corner connection: each of the pair of edges will share ONE zone
    #So e.g. looking at the cutoff = 0.25 / edge city map:
    #154,155 and 154,170 share 154
    thisedge <- get.edgelist(frontier.graph)[i,]
    
    zones <- c(
      as.numeric(unlist(strsplit(thisedge[1],","))),
      as.numeric(unlist(strsplit(thisedge[2],",")))
    )
    
    #If it's NOT a corner connection...
    if(length(unique(zones))==length(zones)){
      
      cat(i," -- ", zones,": straight connection\n")
      
      #We're just going to check if each pair has the same polarity in the same spatial direction
      #So avoiding checkerboard patterns / only keeping connected frontiers with the "same side" 
      
      #Job 1: get pairs on the same side so we know which way polarity should go
      #This will be e.g. one zone in one edge: zone it has rook contiguity with that is
      #(a) NOT in its own pair and (b) IS in the other pair
      #Can get the second "side" directly from that too (! it)
      
      #So pick a random pair.
      randompair <- as.numeric(unlist(strsplit(thisedge[1],",")))
      otherpair <- as.numeric(unlist(strsplit(thisedge[2],",")))
      
      #Find rook contig for that pair that does NOT include its pair but DOES include one from the other pair
      #samesideneighbour and othersideneighbour will be the same numbers as otherpair, but now we'll know WHERE they are
      #In rel to each side
      samesideneighbour <- 
        neigbs.rook[[randompair[1]]][ neigbs.rook[[randompair[1]]] %in% otherpair & !neigbs.rook[[randompair[1]]] %in% randompair  ]
      
      #This is just NOT randompair[1] or its pair, and NOT samesideneighbour
      othersideneighbour <- zones[!zones %in% c(randompair,samesideneighbour)]
      
      cat("same side pairs found: ", randompair[1],",",samesideneighbour, " and ",randompair[2],",",othersideneighbour,"\n",sep = "")
      
      #Now: just find difference between the two pairs but keeping the order spatially correct. 
      #Then multiply – if pos, they're same-side, if neg, they're not.
      pair1diff = sourcemap$valuecol[samesideneighbour] - sourcemap$valuecol[othersideneighbour]
      pair2diff = sourcemap$valuecol[randompair[1]] - sourcemap$valuecol[randompair[2]]
      
      if(pair1diff * pair2diff > 0){
        
        cat('valid straight same-side connection \n')
        
        #Then flag this edge/frontier connection pair as a valid connection
        E(frontier.graph)[i]$validedge <- T
        
      } else {
        
        cat('NOT valid straight same-side connection \n')
        
        #Then flag this edge/frontier connection pair as a valid connection
        E(frontier.graph)[i]$validedge <- F
        
      }
      
      
    } else {
      
      #Else, it's a corner connection
      cat(i," -- ", zones,": corner connection\n")
      
      #Find corner zone
      corner <- as.numeric(names(which.max(table(zones)))[1])
      #And two non-corner zones
      noncorners <- zones[!zones %in% corner]
      
      #We need the corner queen contig zone that is ALSO rook contig to BOTH the other two non-corner zones
      queen.corner <- neigbs.queen[[corner]][neigbs.queen[[corner]] %in% neigbs.rook[[ noncorners[1] ]] ]
      queen.corner <- queen.corner[ queen.corner %in% neigbs.rook[[ noncorners[2] ]] ]
      
      #OK, now we have our four zones.
      #Now we want to make sure all the non-corner zones are the *same* polarity difference
      #Compared to the corner AND above the frontier threshold.
      #We already know the rook contig ones are above the threshold - 
      #can we use that just to then check the queen-corner?
      
      #Only need one of them to test against, actually
      
      #Nah, I think we need both, otherwise polarity can still be wrong
      
      corner_v_noncorner1 = sourcemap$valuecol[corner] - sourcemap$valuecol[noncorners[1]]
      corner_v_noncorner2 = sourcemap$valuecol[corner] - sourcemap$valuecol[noncorners[2]]
      
      corner_v_queencorner = sourcemap$valuecol[corner] - sourcemap$valuecol[queen.corner]
      
      #If the difference is the same polarity, multiplying will give us a positive number
      #ACD also needs to be beyond cutoff
      if(corner_v_noncorner1 * corner_v_queencorner > 0 & 
         corner_v_noncorner2 * corner_v_queencorner > 0 &
         abs(corner_v_queencorner) > frontiercutoff){
        
        cat('valid corner\n')
        
        #Then flag this edge/frontier connection pair as a valid connection
        E(frontier.graph)[i]$validedge <- T
        
      } else {
        
        cat('NOT valid corner\n')
        E(frontier.graph)[i]$validedge <- F
        
      }
      
    }#end corner connection test
    
    
  }#End check connected frontier validity
  
  
  
  #With this current map, those should all be valid... tick
  # table(E(frontier.graph)$validedge)
  
  #Make new graph filtering on the trues
  #Somehow could also do with making a new map from those pairs.
  #Note, we're DELETING EDGES WE DON'T WANT, hence setting to false
  frontier.graph.valid <- delete_edges(frontier.graph, E(frontier.graph)[E(frontier.graph)$validedge==F])
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #ADD NEW VALID FRONTIERS TO FRONTIER SF----
  
  #Now let's try and plot the new valid frontiers
  #NOTE: e.g. in the checkerboard, we're still KEEPING single frontiers
  #But now they're NOT CONNECTED
  #So - to visualise the change, show frontiers in components larger than 1
  
  #So we need some way to add the component size to the vertex (one particular frontier) as an attribute
  membership <- as.data.frame(components(frontier.graph.valid)$membership) %>% 
    rename(membernames = `components(frontier.graph.valid)$membership`)
  
  membership$vertex <- rownames(membership)
  
  componentsize <- as.data.frame(components(frontier.graph.valid)$csize) %>% 
    rename(csize = `components(frontier.graph.valid)$csize`)
  
  componentsize$membernames <- as.numeric(rownames(componentsize))
  
  #Now we can merge those two
  componentsize_n_frontiers <- membership %>% left_join(componentsize, by = 'membernames')
  
    #Actually, we should be able to match against id and id1 in frontiers (found for above map)
  #They're both in the correct order
  
  #Can match directly if I make a new field in frontiers, then.
  frontiers <- frontiers %>% 
    mutate(vertex = paste0(id,",",id.1))
  
  #Merge component size onto the vertices
  frontiers <- frontiers %>%
    left_join(componentsize_n_frontiers, by = 'vertex')
  
  frontiers$csize_factor <- factor(frontiers$csize)
  
  #DONE. 
  #RETURN BOTH FULL FRONTIER GRAPH AND 'VALID' FRONTIER GRAPH
  #PLUS FRONTIER SF
  return(list(frontier.graph = frontier.graph, frontier.graph.valid = frontier.graph.valid, frontiers.sf = frontiers))
  
}




