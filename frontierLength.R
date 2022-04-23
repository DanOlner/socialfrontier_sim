#MEASURE FRONTIER LENGTH
#Meaning here: find out which frontiers (edges between zones with some cutoff for demog difference)
#touch each other
#Make graph of that
#Count number of edges in each subgraph - where each subgraph is one frontier
library(igraph)
library(tidyverse)
library(raster)
library(sf)
library(spdep)
library(tmap)

source('FUNCTIONS_optimiseAACD.R')
Rcpp::sourceCpp("simfunctions.cpp")


#Let's use the two imaginary cities
cities <- read_csv('larger_simcities_forcsv_16x16_3citiesstacked.csv',col_names = F)

#Darn it, goes down rows then across. Of course it does. Harrumph!
#Wouldn't matter if we hadn't joined three together, mind. Still...
cities.c <- as.numeric(unlist(cities))

#More faff required...
smooth.c <- cities[1:16,] %>% as.matrix()
smooth.c <- t(smooth.c)
smooth.c <- as.numeric(smooth.c)

edge.c <- cities[17:32,] %>% as.matrix()
edge.c <- t(edge.c)
edge.c <- as.numeric(edge.c)

edgesplit.c <- cities[33:48,] %>% as.matrix()
edgesplit.c <- t(edgesplit.c)
edgesplit.c <- as.numeric(edgesplit.c)

#Order seems incredibly random in QGIS grid version! How odd.
#Let's try sf instead. Use same shapefile as base to make grid on...
base <- st_read('saves/shapefiles/16by16gridSquare.shp')

#Create new grid over top... ah ha, is now in right order!!
gridsmooth <- st_make_grid(base,n=c(16,16)) %>% st_sf()
gridsmooth$id <- c(1:256)
plot(gridsmooth)

#Repeat for other two. Might want to put in list to process better once made.
gridedge_centre <- st_make_grid(base,n=c(16,16)) %>% st_sf()
gridedge_centre$id <- c(1:256)
plot(gridedge_centre)

gridedge_split <- st_make_grid(base,n=c(16,16)) %>% st_sf()
gridedge_split$id <- c(1:256)
plot(gridedge_split)

#Add data from cities
gridsmooth$peeps1 <- smooth.c
gridedge_centre$peeps1 <- edge.c
gridedge_split$peeps1 <- edgesplit.c

#Check those all contain the same values... tick
table(gridsmooth$peeps1[order(gridsmooth$peeps1)]==gridedge_centre$peeps1[order(gridedge_centre$peeps1)])
table(gridedge_split$peeps1[order(gridedge_split$peeps1)]==gridedge_centre$peeps1[order(gridedge_centre$peeps1)])

#..
#Neighbour lists
#convert to spatialpolygondataframe so we can use spdep
# gridsmooth.spdep <- as_Spatial(gridsmooth)
# gridedge_centre.spdep <- as_Spatial(gridedge_centre)
# gridedge_split.spdep <- as_Spatial(gridedge_split)
# 
# #https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
# #GIVES US A NEIGHBOUR LIST FOR EVERY WARD
# gridsmooth.neighbours <- poly2nb(gridsmooth.spdep)
# gridedge_centre.neighbours <- poly2nb(gridedge_centre.spdep)
# gridedge_split.neighbours <- poly2nb(gridedge_split.spdep)
# 
# 
# gridsmooth.matrix <- nb2mat(gridsmooth.neighbours,style = 'B')
# 
# #Keep only one-way connections so we're not double-processing zone pairs
# #https://stackoverflow.com/questions/26377199/convert-a-matrix-in-r-into-a-upper-triangular-lower-triangular-matrix-with-those
# gridsmooth.matrix[lower.tri(gridsmooth.matrix)] <- 0
# 
# #Back to neighbour list from that?
# #2nd list item contains the neighbours...? Tick.
# x <- mat2listw(gridsmooth.matrix)
# 
# y <- x[[2]]
# y[[1]]
# y[[17]]#Tick, doesn't include 1.


#function that all up
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

# gridsmooth.neighbours <- uppertri.neighbourlist(gridsmooth)
# gridedge_centre.neighbours <- uppertri.neighbourlist(gridedge_centre)
# gridedge_split.neighbours <- uppertri.neighbourlist(gridedge_split)

#Hah, of course those are all the same neighbours! But could be different so worth doing for now.
#And worth keeping difference clear.
# gridsmooth.neighbours[[1]]
# gridsmooth.neighbours[[17]]

#Actually, let's just use one neighbour list. That's less confusing.


#~~~~~~~~~~~~~~~~~~~~
#TESTING WITH SMALLER 4*4 GRID----
#~~~~~~~~~~~~~~~~~~~~

x <- st_make_grid(base,n=c(4,4)) %>% st_sf()
x$id <- c(1:16)
x$value <- runif(16)

tm_shape(x) +
  tm_fill(col = 'value',legend.show = F) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id')


x.spdep <- as_Spatial(x)
neigbs.rook <- poly2nb(x.spdep,queen = F)
neigbs.queen <- poly2nb(x.spdep,queen = T)

neigbs.rook[[1]]
neigbs.queen[[1]]

#Let's test getting all relevant connected frontiers, assuming (5,6) has a frontier
#What queen contig neighbours do 5 and 6 share?
#1,2,9,10. Yeah, they're all the potential linked frontier zones.
bothqueen <- neigbs.queen[[6]][neigbs.queen[[6]] %in% neigbs.queen[[5]]]

#Now, for all of those, we want pair zones for any that have rook contig WITH EACH OTHER 
#(including the original 5 and 6, though not the 5,6 pair itself)
bothqueen <- c(bothqueen,5,6)

#Cycle through each, extract pairs that have rook contig WITH ANY OTHER IN THE SAME LIST
potentiallink_pairlist <- list()
for(i in bothqueen){
  
  rookoverlap <- bothqueen[bothqueen %in% neigbs.rook[[i]]]
  #Exclude the original ones we're checking
  rookoverlap <- rookoverlap[!rookoverlap %in% c(5,6)]
  
  #Now we can add any we found to the potentiallink_pairlist
  #In the actual code we'll check if there's a frontier and add if so
  for(j in rookoverlap){
    
    #Keep same order so duplicates can be removed
    if(j > i){
    potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(i,j)
    } else {
    potentiallink_pairlist[[length(potentiallink_pairlist)+1]] <- c(j,i)
    }
      # ifelse(j > i, c(i,j),c(i,j)) #Not sure why that didn't work
      
      
    
  }
  
}

#Yup, that's our six potential connected frontiers. Bonza.
#Now just to repeat with a filter for which actually have frontiers.
#And do for all pairs.
potentiallink_pairlist <- unique(potentiallink_pairlist)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GET MAPS OF FRONTIERS FOR 16*16 TEST CITIES----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Cos we need to see where they are!
#Cycle through pairs, check ACD cutouff. If above, intersect zones to get line
# plot(density(gridedge_centre$peeps1))

frontiercutoff <- 0.25

#Check just unique rook pairs (doesn't matter which 16*16 grid, all same neighbours)
uniquerookneighbours <- neighbourlist(gridedge_centre,queen = F,uniquepairs = T)

frontiers <- list()

#Can't use that, need the original zone number
# for(neigh in uniquerookneighbours){
#This instead:
#For each zone...
#Have to do minus here if using unique pairs
#as the final bottom-right zone has no neighbours after others took those pairs
for(i in 1:(nrow(gridedge_centre)-1)){
  
  #And for each neighbour of that zone...
  for(neigh in uniquerookneighbours[[i]]){
  
    #If absolute contig difference is above value for defining as frontier....
    if(abs(gridedge_centre$peeps1[i] - gridedge_centre$peeps1[neigh]) > frontiercutoff){
      
      frontiers[[length(frontiers)+1]] <- st_intersection(gridedge_centre[i,],gridedge_centre[neigh,])
      
    }
      
  }
  
}

#One df:
# frontiers <- bind_rows(frontiers)
#https://geocompr.github.io/presentations/attr.html#6 - weird bind_rows error
frontiers <- do.call(rbind, frontiers)

#Not sure those edge zones should be there?
#plot(st_geometry(frontiers))

tm_shape(gridedge_centre) +
  tm_fill(col = 'peeps1',legend.show = F, n=20) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id') +
tm_shape(frontiers) +
  tm_lines(lwd = 5)
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FIND FULL FRONTIER CONNECTION LIST ON 16*16 WITH FRONTIER CUTOFF----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Assemble an edgelist that igraph can turn into a graph with graph_from_edgelist
#Each edge is a frontier
#Specified by the two zone numbers on each side e.g. (5,6)
#Keep same zone number order so each zonepair is a unique edge

#Just use one neighbour list for each of the three same-dimension grids
neighbs.spdep <- as_Spatial(gridedge_centre)
neigbs.rook <- poly2nb(neighbs.spdep,queen = F)
neigbs.queen <- poly2nb(neighbs.spdep,queen = T)

neigbs.rook[[104]]
neigbs.queen[[104]]

#Let's test getting all relevant connected frontiers, assuming (5,6) has a frontier
#What queen contig neighbours do 5 and 6 share?
#1,2,9,10. Yeah, they're all the potential linked frontier zones.

#Using frontier cutoff from prev. map section
# frontiercutoff <- 0.25

#Not making generic as regards the data column for frontier-finding yet
connectedfrontiers_edgelist = list()
#Make sure to empty this before starting
# potentiallink_pairlist <- list()

#Loop 1: across each main zone...
for(mainzone in 1:nrow(gridedge_centre)){
# for(mainzone in 88){
  
  #And each zone pair made from that zone's rook neighbours
  for(mainzone_neighbour in neigbs.rook[[mainzone]]){
  # for(mainzone_neighbour in 104){
  
    cat('mainzone: ',mainzone, ", mainzone_neighbour: ", mainzone_neighbour, "\n")
    
    #Is that top-level pair a frontier (based on our def)?
    #If not, can skip this check for connected frontiers
    if(abs(gridedge_centre$peeps1[mainzone] - gridedge_centre$peeps1[mainzone_neighbour]) > frontiercutoff){
      
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
          if(abs(gridedge_centre$peeps1[j] - gridedge_centre$peeps1[i]) > frontiercutoff){
          
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


#~~~~~~~~~~~~~~~~~~~~~~~~
#Graph from edge list, find average frontier length----
#~~~~~~~~~~~~~~~~~~~~~~~~

edges.as.matrix <- matrix(unlist(connectedfrontiers_edgelist), ncol = 2, byrow = T)
frontier.graph <- graph_from_edgelist(edges.as.matrix, directed=F)

# plot(frontier.graph)
# components(frontier.graph)
# 
# #List - easy to get average frontier length from this
# mean(components(frontier.graph)$csize)



#~~~~~~~~~~~~~~~~~~~~~~~~~~
#REPEAT FOR RANGE OF CUTOFFS, GET AV FRONTIER LENGTH FOR FULL 0-1 RANGE----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

#Should probably function it up really.
#frontier.graph <- frontiernetworks(frontiercutoff = 0.25, spatial.sf = gridedge_centre)

#Good. Now just need to iterate over frontiercutoff

# results <- lapply(
#   seq(from=0.001, to =1, by = 0.02),function(x){
#     frontier.graph <- frontiernetworks(frontiercutoff = x, spatial.sf = gridedge_centre)
#     return(list(cutoff=x, frontiercount = components(frontier.graph)$no, meanfrontierlength = mean(components(frontier.graph)$csize)))
#   }
# )
# 
# final <- matrix(unlist(results), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)
# 
# ggplot(final,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
#   geom_point() +
#   guides(colour = guide_legend(title="frontier count"))


#As suspected, not linear



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#STAGE 2: CONSISTENT NEIGHBOURING FRONTIERS----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
    pair1diff = gridedge_centre$peeps1[samesideneighbour] - gridedge_centre$peeps1[othersideneighbour]
    pair2diff = gridedge_centre$peeps1[randompair[1]] - gridedge_centre$peeps1[randompair[2]]
    
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
    
    corner_v_noncorner1 = gridedge_centre$peeps1[corner] - gridedge_centre$peeps1[noncorners[1]]
    corner_v_noncorner2 = gridedge_centre$peeps1[corner] - gridedge_centre$peeps1[noncorners[2]]
    
    corner_v_queencorner = gridedge_centre$peeps1[corner] - gridedge_centre$peeps1[queen.corner]
    
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

#Tick
# frontier.graph
# components(frontier.graph)
# frontier.graph.valid
# components(frontier.graph.valid)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#VISUALISE NEW VALID FRONTIERS----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
# V(frontier.graph.valid)$name

#Can match directly if I make a new field in frontiers, then.
frontiers <- frontiers %>% 
  mutate(vertex = paste0(id,",",id.1))

#Merge component size onto the vertices
frontiers <- frontiers %>%
  left_join(componentsize_n_frontiers, by = 'vertex')

frontiers$csize_factor <- factor(frontiers$csize)

#Now we can plot those how we like....
tm_shape(gridedge_centre) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id') +
  tm_shape(frontiers) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(frontiers %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'csize_factor',palette='Set1')
  # tm_lines(lwd = 5, col = 'black')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TEST DIFFERENT 16X16 MAPS----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Using functioned-up version of the code above. E.g.
# x <- find.frontier.lengths(frontiercutoff = 0.25,sourcemap = gridedge_centre, valuecol = 'peeps1')
#Which returns a list with the two graphs (all frontiers and 'valid' ones) and the frontier SF

#Note of maps we have already:
#gridsmooth
#gridedge_centre
#gridedge_split

#They all have the same zones, just re-arranged - and so same index of dissimilarity
dissimilarityindex(pop1 = gridsmooth$peeps1, pop2 = 1 - gridsmooth$peeps1)
dissimilarityindex(pop1 = gridedge_centre$peeps1, pop2 = 1 - gridedge_centre$peeps1)
dissimilarityindex(pop1 = gridedge_split$peeps1, pop2 = 1 - gridedge_split$peeps1)

#We want two, maybe three more, optimised:
#chessboard
#smooth from random

#Use one of the existing as the base
#Secondpop isn't used I don't think...
chess <- gridedge_centre
x <- optimiseAverageAbsoluteContiguousDifference(attribute = chess$peeps1, secondpop = 1-chess$peeps1,
                                                            ncol = 16, nrow = 16, maximise = T, cutoff = 0, breakval = 1000000)

chess$peeps1 <- x$attribute
tm_shape(chess) + tm_fill(col = 'peeps1',legend.show = F, n=20)



#A purely random one to get a good mix
#(Which is closer to checkerboard...)
grid.random <- gridedge_centre
grid.random$peeps1 <- sample(x = grid.random$peeps1, size = nrow(grid.random), replace = F)
#Check DI didn't change... tick
dissimilarityindex(pop1 = grid.random$peeps1, pop2 = 1 - grid.random$peeps1)

tm_shape(grid.random) + tm_fill(col = 'peeps1',legend.show = F, n=20)




#And a smooth one from randomised
smooth.from.rnd <- gridedge_centre
smooth.from.rnd$peeps1 <- sample(x = smooth.from.rnd$peeps1, size = nrow(smooth.from.rnd), replace = F)
#Check DI didn't change... tick
dissimilarityindex(pop1 = smooth.from.rnd$peeps1, pop2 = 1 - smooth.from.rnd$peeps1)

x <- optimiseAverageAbsoluteContiguousDifference(attribute = smooth.from.rnd$peeps1, secondpop = 1-smooth.from.rnd$peeps1,
                                                 ncol = 16, nrow = 16, maximise = F, cutoff = 0, breakval = 1000000)

smooth.from.rnd$peeps1 <- x$attribute
tm_shape(smooth.from.rnd) + tm_fill(col = 'peeps1',legend.show = F, n=20)


getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = smooth.from.rnd$peeps1,ncol = 16,nrow = 16,cutoff = 0)


#~~~~~~~~~~~~~
#OK, we got a lot of maps. Now what are we after?
#Let's check chessboard first - that's the one where 'valid' frontiers should mostly all go away

#Keep frontier cutoff for all our tests here
frontiercutoff = 0.98

chess.result <- find.frontier.lengths(frontiercutoff = frontiercutoff, sourcemap = chess, valuecol = 'peeps1')

chess.result$frontier.graph
chess.result$frontier.graph.valid

components(chess.result$frontier.graph)
components(chess.result$frontier.graph.valid)
max(components(chess.result$frontier.graph.valid)$csize)

hist(components(chess.result$frontier.graph)$csize,breaks=20)
hist(components(chess.result$frontier.graph.valid)$csize, breaks=20)

#Average frontier length for both:
mean(components(chess.result$frontier.graph)$csize)
mean(components(chess.result$frontier.graph.valid)$csize)



table(chess.result$frontiers.sf$csize)

#Aaaand
tm_shape(chess) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id') +
  tm_shape(chess.result$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(chess.result$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'black')
  # tm_lines(lwd = 5, col = 'csize_factor',palette='Set1')


#A better way to show this might be...
tm_shape(chess) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(chess.result$frontiers.sf%>% filter(csize > 2)) +
  # tm_shape(chess.result$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'black') +
  tm_shape(chess.result$frontiers.sf %>% filter(csize == 1)) +
  tm_lines(lwd = 5, col = 'red') +
  tm_shape(chess.result$frontiers.sf %>% filter(csize == 2)) +
  tm_lines(lwd = 5, col = 'green') 

#What is largest one?
tm_shape(chess) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(chess.result$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(chess.result$frontiers.sf %>% filter(csize == max(components(chess.result$frontier.graph.valid)$csize))) +
  tm_lines(lwd = 5, col = 'black')



smooth.rnd.result <- find.frontier.lengths(frontiercutoff = frontiercutoff, sourcemap = smooth.from.rnd, valuecol = 'peeps1')

smooth.rnd.result$frontier.graph
smooth.rnd.result$frontier.graph.valid

components(smooth.rnd.result$frontier.graph)
components(smooth.rnd.result$frontier.graph.valid)

mean(components(smooth.rnd.result$frontier.graph)$csize)
mean(components(smooth.rnd.result$frontier.graph.valid)$csize)


table(smooth.rnd.result$frontiers.sf$csize)

tm_shape(smooth.from.rnd) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'peeps1') +
  # tm_text(text = 'id') +
  tm_shape(smooth.rnd.result$frontiers.sf) +
  tm_lines(lwd = 8, col = 'black') +
  tm_shape(smooth.rnd.result$frontiers.sf %>% filter(csize > 1)) +
  # tm_lines(lwd = 5, col = 'black')
  tm_lines(lwd = 4, col = 'csize_factor',palette='Set1')  +
  tm_layout(legend.outside = T)




grid.random.result <- find.frontier.lengths(frontiercutoff = frontiercutoff, sourcemap = grid.random, valuecol = 'peeps1')

grid.random.result$frontier.graph
grid.random.result$frontier.graph.valid

components(grid.random.result$frontier.graph)
components(grid.random.result$frontier.graph.valid)

#Average frontier length for both:
mean(components(grid.random.result$frontier.graph)$csize)
mean(components(grid.random.result$frontier.graph.valid)$csize)

hist(components(grid.random.result$frontier.graph)$csize,breaks=20)
hist(components(grid.random.result$frontier.graph.valid)$csize, breaks=20)


table(grid.random.result$frontiers.sf$csize)


tm_shape(grid.random) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(chess.result$frontiers.sf%>% filter(csize > 2)) +
  # tm_shape(chess.result$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'black') +
  tm_shape(grid.random.result$frontiers.sf %>% filter(csize == 1)) +
  tm_lines(lwd = 5, col = 'red') +
  tm_shape(grid.random.result$frontiers.sf %>% filter(csize == 2)) +
  tm_lines(lwd = 5, col = 'green') 


tm_shape(grid.random) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id') +
  tm_shape(grid.random.result$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(grid.random.result$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'black')


hist(components(grid.random.result$frontier.graph)$csize,breaks=20)
hist(components(grid.random.result$frontier.graph.valid)$csize, breaks=20)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#REPEAT FULL 0-1 RANGE FOR THRESHOLD WITH NEW METHOD----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Just on chessboard for now.
results.chess <- lapply(
  # seq(from=0.05, to =0.95, by = 0.1),function(x){
    seq(from=0.02, to =0.98, by = 0.02),function(x){
    result <- find.frontier.lengths(frontiercutoff = x, sourcemap = chess, valuecol = 'peeps1')
    return(list(cutoff=x, 
                frontiercount = components(result$frontier.graph.valid)$no, 
                meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
  }
)

final.chess <- matrix(unlist(results.chess), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)

final.chess %>% 
  filter(cutoff > 0.1) %>% 
  summarise(mean(mean_frontierlength))

ggplot(final.chess,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
  geom_point() +
  guides(colour = guide_legend(title="frontier count"))

#Pick one example
x <- find.frontier.lengths(frontiercutoff = 0.02, sourcemap = chess, valuecol = 'peeps1')

tm_shape(chess) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  tm_text(text = 'id') +
  tm_shape(x$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(x$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'csize_factor',palette='Set1') +
  tm_layout(legend.outside=T)

table(components(x$frontier.graph.valid)$csize)
mean(components(x$frontier.graph.valid)$csize)


#Does edge-city-centre look so smooth? Nope.
results.edgecentre <- lapply(
  seq(from=0.02, to =0.98, by = 0.02),function(x){
    result <- find.frontier.lengths(frontiercutoff = x, sourcemap = gridedge_centre, valuecol = 'peeps1')
    return(list(cutoff=x, 
                frontiercount = components(result$frontier.graph.valid)$no, 
                meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
  }
)

final.edgecentre <- matrix(unlist(results.edgecentre), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)

final.edgecentre %>% 
  filter(cutoff > 0.1) %>% 
  summarise(mean(mean_frontierlength))

ggplot(final.edgecentre,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
  geom_point() +
  guides(colour = guide_legend(title="frontier count"))

#A bit arbitrary whether those lower values are one or two connected components isn't it?
x <- find.frontier.lengths(frontiercutoff = 0.02, sourcemap = gridedge_centre, valuecol = 'peeps1')

tm_shape(gridedge_centre) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(x$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(x$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'csize_factor',palette='Set1') +
  tm_layout(legend.outside=T)

table(components(x$frontier.graph.valid)$csize)
mean(components(x$frontier.graph.valid)$csize)

#Can I just check one that I think should have been broken?
#170,171 -- 155,171: should have been broken. Was it? Come back to that...
#Ah, can I see on plot? Noooo
plot(x$frontier.graph.valid)

as_edgelist(x$frontier.graph.valid,names = T)


#SMOOOOTH
results.smooth <- lapply(
  seq(from=0.02, to =0.98, by = 0.02),function(x){
    result <- find.frontier.lengths(frontiercutoff = x, sourcemap = smooth.from.rnd, valuecol = 'peeps1')
    return(list(cutoff=x, 
                frontiercount = components(result$frontier.graph.valid)$no, 
                meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
  }
)

final.smooth <- matrix(unlist(results.smooth), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)

final.smooth %>% 
  filter(cutoff > 0.1) %>% 
  summarise(mean(mean_frontierlength, na.rm=T))

ggplot(final.smooth,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
  geom_point() +
  guides(colour = guide_legend(title="frontier count"))


#What does smooth look like for those high early values?
x <- find.frontier.lengths(frontiercutoff = 0.25, sourcemap = smooth.from.rnd, valuecol = 'peeps1')

tm_shape(smooth.from.rnd) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(x$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red',lty=3) +
  tm_shape(x$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'csize_factor',palette='Set1')+
  tm_layout(legend.outside=T)




#SMOOOOTH ORIGINAL MANUAL
results.smooth.manual <- lapply(
  seq(from=0.02, to =0.98, by = 0.02),function(x){
    result <- find.frontier.lengths(frontiercutoff = x, sourcemap = gridsmooth, valuecol = 'peeps1')
    return(list(cutoff=x, 
                frontiercount = components(result$frontier.graph.valid)$no, 
                meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
  }
)

final.smooth.manual <- matrix(unlist(results.smooth.manual), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)

final.smooth.manual %>% 
  filter(cutoff > 0.1) %>% 
  summarise(mean(mean_frontierlength, na.rm=T))

ggplot(final.smooth.manual,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
  geom_point() +
  guides(colour = guide_legend(title="frontier count"))

#Ah yeah, somewhat obviously, 16 is av length




#How goes it with random?
#Changed to deal with NULL return if no frontiers
#Would be nicer to return empty graphs maybe
grid.random <- gridedge_centre
grid.random$peeps1 <- sample(x = grid.random$peeps1, size = nrow(grid.random), replace = F)


results <- lapply(
  seq(from=0.02, to =0.98, by = 0.02),function(x){
    result <- find.frontier.lengths(frontiercutoff = x, sourcemap = grid.random, valuecol = 'peeps1')
    return(list(cutoff=x, 
                frontiercount = components(result$frontier.graph.valid)$no, 
                meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
  }
)

final <- matrix(unlist(results), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)

final %>% 
  filter(cutoff > 0.1) %>% 
  summarise(mean(mean_frontierlength, na.rm=T))


ggplot(final,aes(x = cutoff, y = mean_frontierlength, colour = factor(frontiercount))) +
  geom_point() +
  guides(colour = guide_legend(title="frontier count"))

#Yup, more or less the same pattern as chessboard, which is kinda what we'd expect
#(Though maybe not the same numbers?)

x <- find.frontier.lengths(frontiercutoff = 0.02, sourcemap = grid.random, valuecol = 'peeps1')

tm_shape(grid.random) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) + 
  # tm_text(text = 'id') +
  tm_shape(x$frontiers.sf) +
  tm_lines(lwd = 3, col = 'red') +
  tm_shape(x$frontiers.sf %>% filter(csize > 1)) +
  tm_lines(lwd = 5, col = 'csize_factor',palette='Set1')+
  tm_layout(legend.outside=T)


table(components(x$frontier.graph.valid)$csize)
mean(components(x$frontier.graph.valid)$csize)



#~~~~~~~~~~~~~~~~~~~~~~~~~~
#FRONTIER LENGTH AV DISTRIBUTION FROM RANDOMISED----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

x <- list()

for(i in 1:5){

  rnd <- gridedge_centre
  rnd$peeps1 <- sample(x = rnd$peeps1, size = nrow(rnd), replace = F)

    results <- lapply(
    seq(from=0.02, to =0.98, by = 0.02),function(x){
      result <- find.frontier.lengths(frontiercutoff = x, sourcemap = rnd, valuecol = 'peeps1')
      return(list(cutoff=x, 
                  frontiercount = components(result$frontier.graph.valid)$no, 
                  meanfrontierlength = mean(components(result$frontier.graph.valid)$csize)))
    }
  )
  
  final <- matrix(unlist(results), ncol = 3, byrow = T) %>% data.frame() %>% rename(cutoff = X1, frontiercount = X2,mean_frontierlength = X3)
  
  x[[length(x)+1]] <- final %>% 
    filter(cutoff > 0.1) %>% 
    summarise(mean(mean_frontierlength, na.rm=T))

}


rez <- unlist(x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
#PLOTS FOR GOOGLE DOCS AMBIENT FRONTIERS ARTICLE----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

#Just the basic maps.
x <- tm_shape(chess) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1) +
  tm_layout(legend.outside=T)
  # tm_text(text = 'peeps1')

tmap_save(x,'outputs/maps_for_googledoc/chess.png',width = 5,height=5)


#AACDs for smooth and grid edge
getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridsmooth$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)
getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridedge_centre$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)
getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = chess$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)

#https://stackoverflow.com/questions/5468280/scale-a-series-between-two-points
library(scales)

rescale(
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridedge_centre$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F),
  to = c(
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridsmooth$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F),
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = chess$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)
  )
)

rescale(3, to=c(0,1))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

range01(c(
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridsmooth$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F),
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridedge_centre$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F),
  getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = chess$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)
))


getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = smooth.from.rnd$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)


#So actually, a simple sorting might even give a lower AACD than opt minimising
#UPDATE: YUP. For grids, sorting gives better minimum than minimising through hillclimb
gridsmooth.fromsorting <- gridedge_centre
gridsmooth.fromsorting$peeps1 <- runif(n = 16*16)
gridsmooth.fromsorting$peeps1 <- sort(gridsmooth.fromsorting$peeps1)

tm_shape(gridsmooth.fromsorting) +
  tm_fill(col = 'peeps1',legend.show = F,n = 20) + 
  tm_borders("gray", lwd = 1)

getAverageAbsoluteContiguousDifference_GRIDONLY(attribute = gridsmooth.fromsorting$peeps1,ncol = 16,nrow = 16,cutoff = 0,torus = F)
  
dissimilarityindex(gridsmooth$peeps1,1-gridsmooth$peeps1)
dissimilarityindex(gridedge_centre$peeps1,1-gridedge_centre$peeps1)

