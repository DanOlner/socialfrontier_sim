library(sf)
library(tidyverse)
library(spdep)
library(raster)
library(tmap)
#library(Rcpp)
#source('SegregationIndicesDuncanLee.R')
Rcpp::sourceCpp("simfunctions.cpp")

#~~~~~~
#DI----
#~~~~~~

grid <- raster(nrow = 30, ncol = 30) %>% rasterToPolygons() %>% as("sf")

#Two random data columns
grid$peeps1 <- runif(n = nrow(grid))
grid$peeps2 <- runif(n = nrow(grid))


#Initial DI using R function version
dissimilarity(grid$peeps1, grid$peeps2)

dissimilarityindex(grid$peeps1,grid$peeps2)

optz <- targetDissimilarityIndex(di_target = 1, threshold = 0.0000001, breakval = 1000000, grid$peeps1,grid$peeps2)

dissimilarityindex(optz$pop1,optz$pop2)


#plot result
gridresult = grid
gridresult$peeps1 = optz$pop1
gridresult$peeps2 = optz$pop2
gridresult$peeps2minuspeeps1 = optz$pop2 - optz$pop1

plot(gridresult)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CHECK ON GETTING ABSOLUTE CONTIG DIFFERENCE----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#CODE FOR CREATING TORUS NEIGHBOUR LIST
#WITH NEIGHBOUR INDEX CORRECTLY REFERENCED IN SPATIAL OBJECT
numcol = 12
numrow = 12

#This is just a neighbour list
#cellz <- cell2nb(nrow = numrow, ncol = numcol, type = "queen", torus = TRUE)
#Rook contiguity
cellz <- cell2nb(nrow = numrow, ncol = numcol, torus = TRUE)

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
spatialz$selection[cellz[[15]]] <- 1

#Tick!
tmap_mode('plot')
qtm(spatialz, fill = 'selection')




#Set some random attribute
spatialz$attribute <- runif(n = nrow(spatialz))

#Nema's way of quickly finding ACD:
x <- spatialz$attribute
n <- length(x)
x.diff.matrix <- outer(x,x,FUN="-")
x.diff.matrix <- abs(x.diff.matrix)

## paired differences amongst geographical neighbours
W0 <- nb2mat(cellz, style = 'B')

x.diff.W0 <- x.diff.matrix[W0 == 1]
x.diff.W0 <- unique(x.diff.W0)#Risky if any other values match, no?





#Testing ACD function
getAverageAbsoluteContiguousDifference(spatialz$attribute, ncol = numcol, nrow = numrow)

#Printed out values in C++ function. Check against spatialz:
check = 91
cellz[[check]]#order is... not consistent? Oh, it's up left right down, even if "up" is actually torus-down. Actually, no, it's just not consistent.
spatialz$attribute[check]
spatialz$attribute[cellz[[check]]]

# i: 90, left side torus. Attribute here: 0.833644, attribute neighbour: 0.193821
# i: 90, right side. Attribute here: 0.833644, attribute neighbour: 0.583322
# i: 90, top side. Attribute here: 0.833644, attribute neighbour: 0.128071
# i: 90, bottom side torus. Attribute here: 0.833644, attribute neighbour: 0.237748

#Check against value in R
#It seems right but needs doubling
sum(x.diff.W0)*2

#Not quite clear why that should be happening. I know Nema's code is removing neighbour duplicates
#But I'm not duplicating in the C++ code.
#Oh well - won't make a difference to optimising to a value anyway.



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#OPTIMISING ACD GIVEN A TARGET DI----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ncol = 30
nrow = 30

#Get a target DI to start with. Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")

#Two random data columns
grid$peeps1 <- runif(n = nrow(grid))
grid$peeps2 <- runif(n = nrow(grid))

optz <- targetDissimilarityIndex(di_target = 0.5, threshold = 0.0000001, breakval = 1000000, grid$peeps1,grid$peeps2)

dissimilarityindex(optz$pop1,optz$pop2)

#plot result
gridresult = grid
gridresult$peeps1 = optz$pop1
gridresult$peeps2 = optz$pop2
gridresult$peeps2minuspeeps1 = optz$pop2 - optz$pop1

#Check order of position... yup, correct top-left ordering, all good.
# gridresult$flag = 0
# gridresult$flag[11] = 1

plot(gridresult)


#Which means we can drop the DI-optimised attribute straight into the ACD calculation.
#Can only use one of the population values to maximise ACD.
#getAverageAbsoluteContiguousDifference(gridresult$peeps1, ncol = ncol, nrow = nrow)


#So now for an optimisation routine.
x <- optimiseAverageAbsoluteContiguousDifference(attribute = gridresult$peeps1, ncol = ncol, nrow = nrow, maximise = F, breakval = 10000000)

#going t osave the 100x100 one, took ages
saveRDS(x,'100x100_0point5DI_minimise.rds')

getAverageAbsoluteContiguousDifference(x, ncol = ncol, nrow = nrow)

gridresult$optimisedACD <- x

plot(gridresult[,c('peeps1','optimisedACD')])

x <- tm_shape(gridresult) +
  tm_fill('optimisedACD', palette = 'viridis', n=20) +
  tm_borders(alpha = 0.5) +
  tm_layout(legend.show = F)

tmap_save(x,'outputs/100x100_DI0point5_minimise.png', width = 2000, height = 2000)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SET PROPORTION, TARGET DI, GET SPREAD OF ACD VALUES----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Will also want DI targeting to maintain proportion. So will require some tweaks while optimising.

ncol = 30
nrow = 30

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")

#We want the 'minority' group to be a set proportion.
#Can then make that proportion spatially even - then target DI value.
proportion = 0.4
totalpop = 100000

#Peeps1 can be minority
#Make spatially even to begin with (so DI = 0, no-one need move to make spatially even)
#Split e.g. 40% of totalpop over each cell grid.
#(I am still not worrying about fractions of people)
grid$peeps1 <- rep((proportion * totalpop)/nrow(grid), nrow(grid))
grid$peeps2 <- rep(((1-proportion) * totalpop)/nrow(grid), nrow(grid))

#Nearly but not quite zero! 
round(dissimilarityindex(grid$peeps1,grid$peeps2),12)

#Check we can target a value
optz <- targetDissimilarityIndex(di_target = 0.5, threshold = 0.000000001, breakval = 1000000, grid$peeps1,grid$peeps2)
optz <- targetDissimilarityIndex(di_target = 1, threshold = 0.000000001, breakval = 1000000, grid$peeps1,grid$peeps2)

dissimilarityindex(optz$pop1,optz$pop2)

#What do the proportions end up being? ~5% drift.
sum(optz$pop1)/(sum(optz$pop1)+sum(optz$pop2))



#Keeping population fixed... 
optz <- targetDissimilarityIndex(di_target = 0.25, threshold = 0.000000001, breakval = 1000000, grid$peeps1, grid$peeps2, keepProportion = T)

dissimilarityindex(optz$pop1,optz$pop2)
#Oh and let's just confirm that halving pop values keeps DI identical... Tick.
#Or in fact dividing by any number, which makes sense.
#It's about relative proportions in each zone, yah.
dissimilarityindex(optz$pop1/5,optz$pop2/8)

#Proportion same, tick. That actually seemed faster... Yeah, it totally is, much. Gets to 1 easily.
sum(optz$pop1)/(sum(optz$pop1)+sum(optz$pop2))

# grid$optpop1 = as.integer(optz$pop1)
# grid$optpop2 = as.integer(optz$pop2)
grid$optpop1 = optz$pop1
grid$optpop2 = optz$pop2

grid$subz = grid$optpop2 - grid$optpop1

# plot(grid)
# plot(grid[,'optpop1'])
# plot(grid[,'subz'])


#Min and max straightforward.
min <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                 secondpop = grid$optpop2,
                                                 ncol = ncol, nrow = nrow, maximise = F, breakval = 1000000)

max <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                 secondpop = grid$optpop2,
                                                 ncol = ncol, nrow = nrow, maximise = T, breakval = 1000000)

getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow)

#Function now also swaps the second pop values 
#Though they're not optimised on
#So that DI stays identical
#Just in case that matters for later or for plotting
# grid$optpop1_acdopt = x$attribute
# grid$optpop2_acdopt = x$secondpop
# dissimilarityindex(grid$optpop1,grid$optpop2)
# dissimilarityindex(grid$optpop1_acdopt,grid$optpop2_acdopt)


x = getRepeatedAACDfromPermutedCells(attribute = grid$optpop1,ncol = ncol, nrow = nrow, numreps = 100000)

ggplot(data.frame(x=x), aes(x = x)) +
  geom_density() +
  geom_vline(xintercept = mean(x), colour = 'green')

#And how does that look against the min and max optimised vals?
#Well... it's visible at least!
ggplot(data.frame(x=x), aes(x = x)) +
  geom_density() +
  geom_vline(xintercept = mean(x), colour = 'green') +
  geom_vline(xintercept = getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow), colour = 'blue') +
  geom_vline(xintercept = getAverageAbsoluteContiguousDifference(max$attribute, ncol = ncol, nrow = nrow), colour = 'blue') 
  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HOW MUCH DOES AACD VARY FOR SAME VALUES OF DI?----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Theory: random targeted DI values can have very different morphologies, presumably
#So what's the range of min/max values?
#This might take a while to run for many but let's just check for a few to see...
#Using same values as last section
reps = 20
rez = data.frame(rep = 1:reps, min = -1, max = -1)

for(i in 1:reps){

  optz <- targetDissimilarityIndex(di_target = 0.25, threshold = 0.00000001, breakval = 1000000, grid$peeps1, grid$peeps2, keepProportion = T)
  
  #Oops, forgot to replace opt values
  
  min <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                     secondpop = grid$optpop2,
                                                     ncol = ncol, nrow = nrow, maximise = F, breakval = 1000000)
  
  max <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                     secondpop = grid$optpop2,
                                                     ncol = ncol, nrow = nrow, maximise = T, breakval = 1000000)

  rez$min[i] = getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow)
  rez$max[i] = getAverageAbsoluteContiguousDifference(max$attribute, ncol = ncol, nrow = nrow)

}


#All together...
allz = data.frame(
  type = c(
    rep("dist",length(x)),
    rep("min",reps),
    rep("max",reps)
  ),
  val = c(x,rez$min,rez$max)
    )

ggplot(allz,aes(x = val, colour = type)) +
  geom_density()

range(rez$min)
range(rez$max)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Get min/max/dist for range of DI vals----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

t = proc.time()

ncol = 30
nrow = 30

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")

#We want the 'minority' group to be a set proportion.
#Can then make that proportion spatially even - then target DI value.
proportion = 0.2
totalpop = 100000

#Peeps1 can be minority
#Make spatially even to begin with (so DI = 0, no-one need move to make spatially even)
#Split e.g. 40% of totalpop over each cell grid.
#(I am still not worrying about fractions of people)
grid$peeps1 <- rep((proportion * totalpop)/nrow(grid), nrow(grid))
grid$peeps2 <- rep(((1-proportion) * totalpop)/nrow(grid), nrow(grid))

results = list()

#for(i in seq(from = 0.05, to = 1, by = 0.05)){
for(i in seq(from = 0.05, to = 1, by = 0.05)){
  
  cat(paste0(i,'\n'))
  
  optz <- targetDissimilarityIndex(di_target = i, threshold = 0.00000001, breakval = 1000000, grid$peeps1, grid$peeps2, keepProportion = T)
  
  #Starts at zero. Pop mix is 40%
  #dissimilarityindex(grid$peeps1, grid$peeps2)
  cat(paste0('Targeted DI result:', dissimilarityindex(optz$pop1,optz$pop2),'\n'))
  
  grid$optpop1 = optz$pop1
  grid$optpop2 = optz$pop2
  
  min <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                     secondpop = grid$optpop2,
                                                     ncol = ncol, nrow = nrow, maximise = F, breakval = 1000000)
  
  max <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
                                                     secondpop = grid$optpop2,
                                                     ncol = ncol, nrow = nrow, maximise = T, breakval = 1000000)
  
  minaacd = getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow)
  maxaacd = getAverageAbsoluteContiguousDifference(max$attribute, ncol = ncol, nrow = nrow)
  
  #And random distribution
  x = getRepeatedAACDfromPermutedCells(attribute = grid$optpop1,ncol = ncol, nrow = nrow, numreps = 25000)
  
  results[[length(results)+1]] = list(di = i, minaacd = minaacd, maxaacd = maxaacd, distribution = x)
  
}

saveRDS(results, 'saves/DI_opt_resultset2_02prop.rds')
results <- readRDS('saves/DI_opt_resultset2_02prop.rds')

proc.time() - t

results[[1]]$di

results.df = data.frame(distribution = results[[1]]$distribution, 
                        di = results[[1]]$di,
                        min = results[[1]]$min,
                        max = results[[1]]$max
                        )

for(i in 2:(length(results))){
  
  results.df = rbind(results.df,
    data.frame(distribution = results[[i]]$distribution, 
                          di = results[[i]]$di,
                          min = results[[i]]$min,
                          max = results[[i]]$max
    )
  )
}



#Start with just distributions.
ggplot(results.df %>% filter(di > 0.1), aes(x = distribution, colour = factor(di))) +
  geom_density() +
  guides(colour = F)

ggsave('outputs/DI_target_distributiondensities.png', width = 5, height = 3.7)

ggplot(results.df, aes(y = distribution, x =fct_rev(factor(di)))) +
  geom_boxplot() +
  geom_point(data = results.df %>% filter(!duplicated(results.df$di)), aes(x = factor(di), y = min), colour = 'green') +
  geom_point(data = results.df %>% filter(!duplicated(results.df$di)), aes(x = factor(di), y = max), colour = 'green') +
  coord_flip()

ggsave('outputs/DI_target_dists_plus_minmaxes.png', width = 5, height = 3.7)




#Get 95% percentiles of the AACD random distribution.
#Plot against the mean. Example:
quantile(results.df$distribution[results.df$di==0.05], c(0.025, 0.975))

mean_n_quantz <- results.df %>% 
  group_by(di) %>% 
  summarise(mean = mean(distribution), `2.5%` = quantile(distribution,0.025), `97.5%` = quantile(distribution,0.975))


#plot that shizzle. It's going to look at bit odd, isn't it?
#Can try two different ways to plot this. 
#1: lines
ggplot(mean_n_quantz %>% gather(key = type, value = value, mean:`97.5%`), aes(x = di, y = value, colour = type)) +
  geom_line() +
  xlab('dissimilarity index') 

ggsave('outputs/AACD_v_DI_lines.png', width = 6, height = 4)

#That's probably the neatest way, but let's look at the other option. Direct plotting:
ggplot(mean_n_quantz, aes(x = di, y = mean)) +
  geom_point() +
  geom_errorbar(aes(x = di, ymin = `2.5%`, ymax = `97.5%`)) +
  xlab('dissimilarity index') +
  ylab('AACD mean & 95% quantiles')

ggsave('outputs/AACD_v_DI_errorbars.png', width = 6, height = 4)



saveRDS(results.df, 'saves/results1.rds')
results2.df <- readRDS('saves/results1.rds')

#Get single df with distribution means
meanz <- results.df %>% 
  group_by(di) %>% 
  summarise(distmean = mean(distribution),
            min = max(min),
            max = max(max))

meanz2 <- results2.df %>% 
  group_by(di) %>% 
  summarise(distmean = mean(distribution),
            min = max(min),
            max = max(max))


meanz$prop = 0.2
meanz2$prop = 0.4

both <- rbind(meanz,meanz2)

#Check
ggplot(meanz %>% 
         gather(key = thing, value = value, distmean:max), 
       aes(y = value, x =fct_rev(factor(di)), colour = thing)) +
  geom_point()

ggplot(both %>% 
         gather(key = thing, value = value, distmean:max), 
       aes(y = value, x =fct_rev(factor(di)), colour = thing)) +
  geom_point() +
  facet_wrap(~prop)

#Check if means of each value are roughly double of each other
both %>% 
  gather(key = thing, value = value, distmean:max) %>% 
  group_by(prop,thing) %>%
  summarise(mean = mean(value))
    
#So linear for proportion. Which is what you'd expect, right?
#Except of course... You can multiply either count of group by any number overall and get the same
#DI but it won't be the same AACD
  
#Tick. Save
write_csv(meanz,'saves/di_vs_aacd.csv')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SMALLER TOY EXAMPLE BASED ON GWILYM'S LITTLE CITY PLOTS----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#From sheet3(3) in Relative Social Frontiers v1.xlsx
#Seven columns eight rows.
ncol = 7
nrow = 8

#Set all this up with the same geography.
#We have two versions to do.
gridsmooth <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")
gridsmooth$id = 1:nrow(gridsmooth)

gridedges <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")
gridedges$id = 1:nrow(gridedges)

#Coding values directly from G (left to right)
#G uses single value - for DI for two pops, the other peeps are 1-x
gridsmooth$peeps1 <- c(0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.6,0.6,0.05,0.2,0.4,0.4,0.4,0.4,0.4,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
gridsmooth$peeps2 <- 1 - c(0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.6,0.6,0.05,0.2,0.4,0.4,0.4,0.4,0.4,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
plot(gridsmooth)

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
plot(gridedges)

#What's current AACD for those two?
smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff=0)
edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow,cutoff=0)

actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD))

#Randomly permute grid cells, get AACD from permutes.
#Smoothperm. Nice. 
smoothperm = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff=0)
edgeperm = getRepeatedAACDfromPermutedCells(attribute = gridedges$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff=0)

#Plot details for those two.
both <- data.frame(smooth = smoothperm, edge = edgeperm) %>% 
  gather(key = type, value = AACD)

mean_n_quantz <- both %>% 
  group_by(type) %>% 
  summarise(mean = mean(AACD), sd = sd(AACD), `2.5%` = quantile(AACD,0.025), `97.5%` = quantile(AACD,0.975))

ggplot(mean_n_quantz, aes(x = type, y = mean)) +
  geom_point() +
  geom_errorbar(aes(x = type, ymin = `2.5%`, ymax = `97.5%`)) +
  xlab('grid arrangement') +
  ylab('AACD permute: mean & 95% quantiles') +
  coord_cartesian(ylim = c(0,1.5)) +
  geom_point(data = actualAACDs, aes(x = type, y = AACD), size = 2, colour = 'blue')

ggsave('outputs/smooth_v_edges_w_actualAACDvalues.png', width = 6, height = 4)


#Actually, that plot - we only need one of the mean/error bars, they're
#both from permuting the same grid squares
#So would be tidier to do something like this:
allthree = data.frame(AACD = c('permute mean','smooth','edge'), 
                      value = c(
                        mean_n_quantz$mean[1]
                      ))



#Get AFPI / equivalent of Z score
#AFPI = [AFI - E(AFI)]/sd(afi)
AFPI <- mean_n_quantz %>% left_join(actualAACDs, by = 'type') %>% 
  mutate(AFPI = abs((AACD - mean)/sd))

#Yup, that looks right.

# min <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
#                                                    secondpop = grid$optpop2,
#                                                    ncol = ncol, nrow = nrow, maximise = F, breakval = 1000000)
# 
# max <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$optpop1, 
#                                                    secondpop = grid$optpop2,
#                                                    ncol = ncol, nrow = nrow, maximise = T, breakval = 1000000)
# 
# minaacd = getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow)
# maxaacd = getAverageAbsoluteContiguousDifference(max$attribute, ncol = ncol, nrow = nrow)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TOY EXAMPLE G's CITY PLOTS: WORKING WITH ACD CUTOFF VALUE----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Amended code in the AACD C++ function so it take a cutoff value
#tricky part: finding the cutoff. 
#That might be easier done using vanilla R code
#So can also plot the maps, see what zones the cutoff value cuts off.




#Actually: just going to find overall diff between min and max as guide to range of ACD min/max
min(gridsmooth$peeps1)
max(gridsmooth$peeps1)
min(gridedges$peeps1)
max(gridedges$peeps1)
max(gridsmooth$peeps1)-min(gridsmooth$peeps1)
max(gridedges$peeps1)-min(gridedges$peeps1)

#Both 0.75 max diff. Min will be zero.


#Now adding cutoff. test effect of cutoff values on AACD overall
#(I'd expect it to rise, right? But needs comparing with permute for same cutoff)
actualresults = list()
permresults = list()

for(cutoff in seq(from=0,to=0.75,by=0.05)){
  
  print(cutoff)

  smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff)
  edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff)
  
  actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD), cutoff = cutoff)
  
  #Randomly permute grid cells, get AACD from permutes.
  #Smoothperm. Nice. 
  smoothperm = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff)
  edgeperm = getRepeatedAACDfromPermutedCells(attribute = gridedges$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff)
  
  #Plot details for those two.
  both <- data.frame(smooth = smoothperm, edge = edgeperm) %>% 
    gather(key = type, value = AACD)
  
  mean_n_quantz <- both %>% 
    group_by(type) %>% 
    summarise(mean = mean(AACD), sd = sd(AACD), `2.5%` = quantile(AACD,0.025), `97.5%` = quantile(AACD,0.975)) %>% 
    mutate(cutoff = cutoff)
  
  actualresults[[length(actualresults)+1]] <- actualAACDs
  permresults[[length(permresults)+1]] <- mean_n_quantz

}

actualallz <- bind_rows(actualresults)
permallz <- bind_rows(permresults)


ggplot(permallz, aes(x = type, y = mean)) +
  geom_point() +
  geom_errorbar(aes(x = type, ymin = `2.5%`, ymax = `97.5%`)) +
  xlab('grid arrangement') +
  ylab('AACD permute: mean & 95% quantiles') +
  coord_cartesian(ylim = c(0,1.5)) +
  geom_point(data = actualallz, aes(x = type, y = AACD), size = 2, colour = 'blue') +
  facet_wrap(~cutoff, nrow = 1)

ggsave('outputs/smooth_v_edges_w_actualAACDvalues_cutoffs.png', width = 12, height = 4)





#~~~~~~~~~~~~~~~~~~~~~~~~
#VISUALISE BORDER CUTOFFS----
#~~~~~~~~~~~~~~~~~~~~~~~~

#To be clear on how the cutoff is changing what gets averaged

#Using the data above, make a sp object.
#gwilym's zones are 8 rows, 7 columns.
spatialz <- raster(nrow = 8, ncol = 7) %>% rasterToPolygons() %>% as("sf")
spatialz$smooth = gridsmooth$peeps1
spatialz$edge = gridedges$peeps1

#Yup, that's it.
plot(spatialz)

#Still using predefined neighbour list that includes torus
#Actually, I should probably stop with torus: That would mess with G's intention with those grids.
# cellz <- cell2nb(nrow = nrow, ncol = ncol, torus = FALSE, type = 'rook')
#  
# #Look, torus!
# cellz[[2]]
# 
# ## paired differences amongst geographical neighbours
# W0 <- nb2mat(cellz, style = 'B')
# 
# #Get paired absolute differences (Nema's code in permutation test)
# x.diff.matrix <- outer(gridsmooth$peeps1,gridsmooth$peeps1,FUN="-")
# x.diff.matrix <- abs(x.diff.matrix)
# 
# #And in theory this should work, though I recall something about the cellz creation being 
# #The wrong orientation?
# x.diff.W0 <- x.diff.matrix[W0 == 1]



#x.diff.W0 <- unique(x.diff.W0)#Risky if any other values match, no? Yup, won't work here!

#Grab meng le's code for getting borders into an sf
#Err, then shouldn't need the above
#But it won't reflect 

#Except the intersection approach won't work for torus edges.
#But maybe I don't need those.
#But then the issue is getting the right intersections.



#Argh. Tempted to think it's going to be easier to do this.
#Note: not doing for torus cos intersection won't work round torus.
borders = list()

for(x in 1:length(cellz)){
  
  for(nb in cellz[[x]]){
    
    border <- st_intersection(
      spatialz[x,],
      spatialz[nb,]
    )
    
    #border$acd <- abs(spatialz$smooth[x] - spatialz$smooth[nb])
    
    borders[[length(borders)+1]] <- border
    
  }
  
}


spatialz2 <- st_transform(spatialz, 27700)

border <- st_intersection(
  spatialz2[2,],
  spatialz2[3,]
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#USE NEIGHBOUR LISTS IN RCPP----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#First I need some actual neighbours
#Quick check, with torus
numcol = 12
numrow = 12

#Rook contiguity
cellz <- cell2nb(nrow = numrow, ncol = numcol, torus = TRUE)

cellz[[1]]

#zero is 1!
testNeighbours(cellz,0)

#Bonza, we can access neighbours for any neighbour list
#I imagine this'll slow things down a bit, but a good start
#Given we need a specific function.
#Although there's prob a matrix math way of applying that function quickly
#RCPParmadillo? Will see.






