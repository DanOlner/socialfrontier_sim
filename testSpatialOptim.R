#Testing optim for maximising/minimising spatial layouts
library(sf)
library(tidyverse)
library(spdep)
library(raster)
library(tmap)
source('functions_r.R')

#~~~~~~
#DI----
#~~~~~~

grid <- raster(nrow = 30, ncol = 30) %>% rasterToPolygons() %>% as("sf")

#Two random data columns
grid$peeps1 <- runif(n = nrow(grid))
grid$peeps2 <- runif(n = nrow(grid))

#Initial DI:
dissimilarity(grid$peeps1, grid$peeps2)

#copy to optimise
gridopt <- grid

#Randomly replace cells with new values. Attempt to increase DI.
#Test first with fixed run.
for(i in 1:30000){
  
  #hill climb copy, dump if no improvement
  gridoptcopy <- gridopt
  
  if(runif(1) > 0.5){
    gridoptcopy$peeps1[sample(1:nrow(gridopt),size = 1)] <- runif(1)
  } else {
    gridoptcopy$peeps2[sample(1:nrow(gridopt),size = 1)] <- runif(1)
  }
  
  if(
    dissimilarity(gridoptcopy$peeps1, gridoptcopy$peeps2) >
    dissimilarity(gridopt$peeps1, gridopt$peeps2)
  ){
    gridopt = gridoptcopy
    # cat('increased: ',dissimilarity(gridoptcopy$peeps1, gridoptcopy$peeps2),'\n')
  } else {
    # cat('no increase\n')
  }
  
  if(i %% 500 == 0) cat(i,'\n')
  
}

#DI after 30K rounds? 0.9472624
dissimilarity(gridopt$peeps1,gridopt$peeps2)

plot(grid[,'peeps1'])
plot(gridopt[,'peeps1'])

saveRDS(grid,'saves/startinggrid_DI.rds')
saveRDS(gridopt,'saves/optimisedgrid_DI.rds')


#~~~~~~~~~~~~~
#MORAN'S I----
#~~~~~~~~~~~~~

#Ooo, what is this spdep grid that can do a torus??
#https://www.rdocumentation.org/packages/spdep/versions/1.0-2/topics/cell2nb

numcol = 30
numrow = 30

#This is just a neighbour list
cellz <- cell2nb(nrow = numrow, ncol = numcol, type = "queen", torus = TRUE)

#Look, torus!
cellz[[1]]

#I'm guessing we can't actually get a spatial object from that?

#via https://cran.r-project.org/web/packages/spdep/spdep.pdf
xyc <- attr(cellz, "region.id")
xy <- matrix(as.integer(unlist(strsplit(xyc, ":"))), ncol=2, byrow=TRUE)
#plot(cellz, xy)

#Think we need to make an equivalent spatial object
spatialz <- raster(nrow = numrow, ncol = numcol) %>% rasterToPolygons() %>% as("sf")

#And just make sure they match...
spatialz$row <- xy[,2]
spatialz$col <- xy[,1]

spatialz$selection <- 0
spatialz$selection[cellz[[320]]] <- 1

#Tick!
tmap_mode('view')
qtm(spatialz, fill = 'selection')



#OPTIMISE MORAN FOR THAT

#list for moran.test
list4opt <- nb2listw(cellz)

#For Duncan's function
W <- nb2mat(cellz)


for(j in 1:5){

  #Start with random values
  spatialz$value <- runif(n = nrow(spatialz))
  spatialzopt <- spatialz

  for(i in 1:10000){
    
    spatialzoptcopy <- spatialzopt
    
    currentmoran <- I.compute(spatialzopt$value,1,W)
    
    spatialzoptcopy$value[sample(1:nrow(spatialzoptcopy),1)] <- runif(1)
    
    if(
      #moran.test(x = spatialzoptcopy$value, listw = list4opt)$estimate[1] >
      #moran.test(x = spatialzopt$value, listw = list4opt)$estimate[1] 
      
      #Duncan Lee's function is faster to just get basic value. Total cell value is 1 in this case.
      I.compute(spatialzoptcopy$value,1,W) > currentmoran
      
    ){
      spatialzopt = spatialzoptcopy
    } 
    
    if(i %% 10 == 0){
      cat(i,': ',currentmoran,'\n')
    }
    
      
  }#end for i
  
  both <- rbind(
    spatialz %>% mutate(state = 'initial'),
    spatialzopt %>% mutate(state = 'optimised')
  )
  
  x <- tm_shape(both) +
    tm_fill('value', palette = 'viridis') +
    tm_borders(alpha = 0.5) +
    tm_layout(legend.show = F) +
    tm_facets(by = 'state')
  
  save_tmap(x, paste0('saves/optimised_morans/',j,'.png'))
  
  
}#end for j


#Is the split about 50/50?
ggplot(spatialzopt, aes(x = value)) +
  geom_density()
ggplot(spatialzopt, aes(x = value)) +
  geom_histogram()


moran.test(x = spatialz$value, listw = list4opt)$estimate[1]
I.compute(spatialz$value,1,W)
moran.test(x = spatialzopt$value, listw = list4opt)$estimate[1]

plot(spatialz[,'value'])
plot(spatialzopt[,'value'])

tmap_mode('plot')

tm_shape(spatialz) +
  tm_fill('value', palette = 'viridis') +
  tm_borders(alpha = 0.5) +
  tm_layout(legend.show = F)

tm_shape(spatialzopt) +
  tm_fill('value', palette = 'viridis') +
  tm_borders(alpha = 0.5) +
  tm_layout(legend.show = F)

both <- rbind(
  spatialz %>% mutate(state = 'initial'),
  spatialzopt %>% mutate(state = 'optimised')
              )

x <- tm_shape(both) +
  tm_fill('value', palette = 'viridis') +
  tm_borders(alpha = 0.5) +
  tm_layout(legend.show = F) +
  tm_facets(by = 'state')

save_tmap(x, 'saves/optimised_morans/x.png')



#speed diff
x <- proc.time()
for(i in 1:20) moran.test(x = spatialz$value, listw = list4opt)$estimate[1]
proc.time() - x

#Yup, much faster
x <- proc.time()
for(i in 1:20) I.compute(spatialz$value,1,W)
proc.time() - x



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MAXIMISE MORAN, MINIMISE DI----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Using 1-value to get two values for dissimilarity

for(j in 1:5){
  
  #Start with random values
  spatialz$value <- runif(n = nrow(spatialz))
  spatialzopt <- spatialz
  
  for(i in 1:20000){
    
    spatialzoptcopy <- spatialzopt
    
    currentmoran <- I.compute(spatialzopt$value,1,W)
    #So optimising these would other things being equal push cell values to 0.5...?
    currentDI <- dissimilarity(spatialzopt$value, 1 - spatialzopt$value)
    
    spatialzoptcopy$value[sample(1:nrow(spatialzoptcopy),1)] <- runif(1)
    
    if(
      #moran.test(x = spatialzoptcopy$value, listw = list4opt)$estimate[1] >
      #moran.test(x = spatialzopt$value, listw = list4opt)$estimate[1] 
      
      #Duncan Lee's function is faster to just get basic value. Total cell value is 1 in this case.
      #Allow to be equal so things can continue of one hits a peak
      I.compute(spatialzoptcopy$value,1,W) >= currentmoran &
      dissimilarity(spatialzoptcopy$value, 1 - spatialzoptcopy$value) <= currentDI
      
    ){
      spatialzopt = spatialzoptcopy
    } 
    
    if(i %% 10 == 0){
      cat(i,': moran = ',currentmoran,'   DI = ',currentDI,'\n')
    }
    
    
  }#end for i
  
  both <- rbind(
    spatialz %>% mutate(state = 'initial'),
    spatialzopt %>% mutate(state = 'optimised')
  )
  
  x <- tm_shape(both) +
    tm_fill('value', palette = 'viridis') +
    tm_borders(alpha = 0.5) +
    tm_layout(legend.show = F) +
    tm_facets(by = 'state')
  
  save_tmap(x, paste0('saves/optimised_moran_n_DI/',j,'.png'))
  
  
}#end for j






