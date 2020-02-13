library(sf)
library(tidyverse)
library(spdep)
library(raster)
library(tmap)
library(stringr)
library(cowplot)
theme_set(theme_grey())
source('FUNCTIONS_optimiseAACD.R')

# library(RcppArmadillo)
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
#100x100 grid: even spread for proportion ethnicity----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#So that there are proportion values spread out evenly from 0 to 1
#So it's possible to get smooth surfaces
#Minimising AACD for this kind of DI should give us that
ncol = 50
nrow = 50

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")


grid$peeps1 <- runif(nrow(grid))
grid$peeps2 <- 1-grid$peeps1
hist(grid$peeps1)#Evenly spread from 0 to 1

#NOISE
plot(grid[,'peeps1'])


#Optimise to minimum
min <- optimiseAverageAbsoluteContiguousDifference(attribute = grid$peeps1, 
                                                   secondpop = grid$peeps2,
                                                   ncol = ncol, nrow = nrow, maximise = F, breakval = 5000000, cutoff=0)



class(min)
grid$minimised <- min$attribute
plot(grid[,'minimised'])

#So that's the pattern we're familiar with
# saveRDS(grid,'saves/grid100x100_uniformprop_w_minimised.rds')


#Put more colours in
x <- tm_shape(grid) +
  tm_fill("minimised", n=30, legend.show = F)

tmap_save(tm = x, filename = 'outputs/50x50_DI_uniform_minimise.png', width=10,height=10)


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

#Check on cutoffs
cutoffz <- gridsmooth
cutoffz$peeps1 <- ifelse(cutoffz$peeps1 < 0.25,-1,cutoffz$peeps1)
plot(cutoffz)

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
plot(gridedges)

#What's current AACD for those two?
Rcpp::sourceCpp("simfunctions.cpp")
# smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff=0,torus=T)
# edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow,cutoff=0,torus=T)

#Check original version is same, new one with torus option didn't break. Tick.
#ORIG_getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff=0)
#ORIG_getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow,cutoff=0)

#Nontorus
smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff=0, torus=F)
edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow,cutoff=0,torus=F)



#While here: check using neighbours finds the same result
##Get neighbour list
gridedges.sp <- as_Spatial(gridedges)
gridsmooth.sp <- as_Spatial(gridsmooth)

#Rook contig plz! We only want bordering cells
gridedges.neighbours <- poly2nb(gridedges.sp, queen=F)
gridsmooth.neighbours <- poly2nb(gridsmooth.sp, queen=F)
# 
# #TICK
x <- getNeighbourIndexAACD(attribute = gridsmooth$peeps1,nblist = gridsmooth.neighbours)
y <- getNeighbourIndexAACD(attribute = gridedges$peeps1,nblist = gridedges.neighbours)


hist( x$allACDs , xlim = c(0,1))
hist(y$allACDs, xlim = c(0,1))


actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD))

#Randomly permute grid cells, get AACD from permutes.
#Smoothperm. Nice.
#mean is roughly...(same for both cos same cells values just repositioned)
#Torus not huge diff to randomised
mean(getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 2000, cutoff=0, torus=F))

smoothperm = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff=0, torus=F)
edgeperm = getRepeatedAACDfromPermutedCells(attribute = gridedges$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff=0, torus=F)

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

ggsave('outputs/smooth_v_edges_w_actualAACDvalues_nontorus.png', width = 6, height = 4)


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


#Try edge detection
#https://www.rdocumentation.org/packages/wvtool/versions/1.0/topics/edge.detect
library(wvtool)

#need it to be a raster, right?
gridsmooth_raster <- raster(nrow = nrow, ncol = ncol)
values(gridsmooth_raster) <- gridsmooth$peeps1
plot(gridsmooth_raster)

gridsmooth_raster <- rgb2gray(gridsmooth_raster)
rez <- edge.detect(gridsmooth_raster, thresh1=1, thresh2=15, noise="gaussian", noise.s=3, method="Canny")




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

  smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
  edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
  
  actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD), cutoff = cutoff)
  
  #Randomly permute grid cells, get AACD from permutes.
  #Smoothperm. Nice. 
  smoothperm = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff,torus=F)
  edgeperm = getRepeatedAACDfromPermutedCells(attribute = gridedges$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff, torus=F)
  
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

ggsave('outputs/smooth_v_edges_w_actualAACDvalues_cutoffs_notorus.png', width = 12, height = 4)





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

Rcpp::sourceCpp("simfunctions.cpp")

#First I need some actual neighbours
#Quick check, with torus
numcol = 12
numrow = 12

#Rook contiguity
cellz <- cell2nb(nrow = numrow, ncol = numcol, torus = TRUE)

sapply(1:length(cellz), function(x) cellz[[x]])

cellz[[10]]

#zero is 1!
testNeighbours(cellz,10)

#Bonza, we can access neighbours for any neighbour list
#I imagine this'll slow things down a bit, but a good start
#Given we need a specific function.
#Although there's prob a matrix math way of applying that function quickly
#RCPParmadillo? Will see.

#Actually, as long as we're sticking to a grid, doing this manually is going to be faster.

#random attribute column for that data
#Just an arbitrary neighbour list, doesn't matter what
attrs <- runif(length(cellz))

weightedNeighbourIndex(attrs,cellz)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TEST CONTIG MATRIX IN ARMADILLO----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rcpp::sourceCpp("armadillo.cpp")

#Get contig matrix from above torus
x <- nb2mat(cellz)#default is row-normalised

#Just try a random vector for multiplying
y = runif(numrow*numcol)

results <- armatest(x,y)

#Is that correct? Tick.
baserez <- x %*% y
table(results == baserez)


#Updated rcpp to loop 1000 times. Check speed diff
repbase <- function(x,y){
  for(i in 1:5000) x %*% y
}

library(rbenchmark)
benchmark(repbase(x,y), armatest(x,y))#3.55' vs 13.22' for 5000 reps. 3.5 times faster than base. Not vast.


#See about using sparse matrix
#Yup, this gets it working
#https://stackoverflow.com/questions/10555210/r-convert-matrix-or-data-frame-to-sparsematrix
A <- as(x, "sparseMatrix")#Both the same?
A <- Matrix(x, sparse = TRUE) 

res_sparse <- armatest2(A,y)#Tick!
table(res_sparse == results)

#Speed? Putting back as 5000 loop with void return
#Sparse faster but difference not huge 1:1.2 ratio, ish.
benchmark(armatest(x,y), armatest2(A,y))





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CHECK ON NEMA'S CODE FOR FINDING ABS CONTIG DIFF VIA MATRIX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#On some small example
numcol = 4
numrow = 4

#Rook contiguity
cellz <- cell2nb(nrow = numrow, ncol = numcol, torus = TRUE)

cellz[[1]]

#zero is 1!
testNeighbours(cellz,0)

x <- nb2mat(cellz)#default is row-normalised

#Just try a random vector for multiplying
y = runif(numrow*numcol)


#Abs contig diffs for all those pairs
x.diff.matrix <- outer(y,y,FUN="-")
x.diff.matrix <- abs(x.diff.matrix)

## paired differences amongst geographical neighbours
# x.diff.W0 <- x.diff.matrix[x == 1]
#Oops - row normalised, so...
x.diff.W0 <- x.diff.matrix[x > 0]

#this part I worry about. What if two are the same?
# x.diff.W0 <- unique(x.diff.W0)


#Feel like halving it might be better/easier.
#Just find mean then halve.
mean(x.diff.W0)/2
mean(unique(x.diff.W0))
#Oh of course: THE MEAN WILL BE THE SAME:
#It's twice the count of cells and each value is repeated twice.
#So don't need to divide by two and actually don't need unique either
mean(x.diff.W0)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#WEIGHTED NEIGHBOUR INDEX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#See notes section "new measure?"
#Use G's little sample cities again - will maybe scale them up at some point
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

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)


#Get neighbour list
gridedges.sp <- as_Spatial(gridedges)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(gridedges.sp, queen=F)

#Tick, max four
neighbours[[25]]

displayAllNeighbours(gridedges$peeps1,neighbours)

#CRASH
#Do same for loops in R, look for issues
#This is the key line:
#all we're looking for: j neighbours of neighbours (2nd order) that are also i neighbours (first order).
#(NOT overlapping neighbours of neighbours of j and neighbours of neighbours of i - they don't overlap!)
for(i in 1:length(gridedges$peeps1)){
  
  #Get neighbours of cell i
  i_nb <- neighbours[[i]]
  
  #For each neigbour of cell i: j
  for(j in 1:length(i_nb)){
    
    #Get neighbours of cell j
    j_nb <- neighbours[[ i_nb[[j]] ]]
    
    #We then want neighbours of neighbours of cell j
    #all we're looking for: j neighbours of neighbours (2nd order) that are also i neighbours (first order).
    j_nb_of_nb <- lapply(j_nb, function(x) neighbours[[x]])
    names(j_nb_of_nb) <- j_nb
    j_nb_of_nb <- unlist(j_nb_of_nb)
    
    #What's in both?
    #If we exclude j itself, I think this works
    i_nb[which(i_nb %in% j_nb_of_nb)]
    
    cat(paste0("i:",i," j:",i_nb[[j]],": "))
    cat(i_nb[which(i_nb %in% j_nb_of_nb)])
    cat("\n")
    
  }
  
}



#Currently returning a vector with one result per border
x <- weightedNeighbourIndex(gridedges$peeps1,neighbours)
Rcpp::sourceCpp("simfunctions.cpp")

#testing new
x <- selectiveNeighbourIndex(gridedges$peeps1,neighbours)


#Tick: the c++ script is returning the right number of borders
length(x)
length(unlist(neighbours))

#Which should make it easy-ish to check it's doing what we think it is.
#First entry will be border between one and two
neighbours[[1]]

#result is 0.1
#Looking at diagram: should be average of -
#ACD 1:2, ACD 8:9
#Tick. So that's halving the ACD we'd usually get (0.2) with another border that's 0.
mean(c(abs(gridedges$peeps1[1]-gridedges$peeps1[2]),abs(gridedges$peeps1[8]-gridedges$peeps1[9])))

#Let's pick another with two pairs either side. 17:10 should do.
#Err. How far through is that? That's a bit fiddly to check. 
#Except: use the neighbour list
neighbours[[17]]

#Where in the unlisted do we get that sequence?
nb_vec <- unlist(neighbours)
which(nb_vec==10)
which(nb_vec==16)
which(nb_vec==18)

#53,54,55 are sequential for those borders. So 53 should be border for 17:10
x[53]

#Base ACD would be 17:10 and the pairs each side
mean(c(abs(gridedges$peeps1[17]-gridedges$peeps1[10]),
       abs(gridedges$peeps1[9]-gridedges$peeps1[16]),
       abs(gridedges$peeps1[11]-gridedges$peeps1[18])))

#Tick. Looking promising.

#So, what do we get for our two simple examples?
x <- weightedNeighbourIndex(gridedges$peeps1,neighbours)
y <- weightedNeighbourIndex(gridsmooth$peeps1,neighbours)

mean(x)
mean(y)

#compare to orig
getNeighbourIndexAACD(gridedges$peeps1,neighbours)['AACD']
getNeighbourIndexAACD(gridsmooth$peeps1,neighbours)['AACD']


#Need to compare that to null dist, obv
nullz <- list()
for(i in 1:10000){
  nullz[[length(nullz)+1]] <- mean(weightedNeighbourIndex(sample(gridedges$peeps1,length(gridedges$peeps1),replace = F),
                                                          neighbours))
}

nullz <- unlist(nullz)

ggplot(data.frame(null = nullz), aes(x=nullz)) +
  geom_density() +
  geom_vline(xintercept = mean(x), colour='red') +
  geom_vline(xintercept = mean(y), colour='green')


#Ooo that looks quite promising!
#Now we want to maximise / minimise to see what we get
#For which we could do with setting up something generic given it's only returning a single value
#Though for now...
Rcpp::sourceCpp("simfunctions.cpp")

#Check... tick
getWeightedNeighbourIndexMean(gridedges$peeps1,neighbours)
mean(x)



#TEST USING JUST ONE IJ PAIR PER CELL
Rcpp::sourceCpp("simfunctions.cpp")
getWeightedNeighbourIndexMeanKeepOnlyMaxACDperSquare(gridedges$peeps1,neighbours)
getWeightedNeighbourIndexMeanKeepOnlyMaxACDperSquare(gridsmooth$peeps1,neighbours)

#And what's random?
getWeightedNeighbourIndexMeanKeepOnlyMaxACDperSquare(sample(gridsmooth$peeps1,length(gridsmooth$peeps1),replace = F),neighbours)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#OPTIMISE WEIGHTED WITH UNIFORM 50*50 GRID----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#So that there are proportion values spread out evenly from 0 to 1
#So it's possible to get smooth surfaces
#Minimising AACD for this kind of DI should give us that
ncol = 20
nrow = 20

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")


grid$peeps1 <- runif(nrow(grid))
grid$peeps2 <- 1-grid$peeps1
hist(grid$peeps1)#Evenly spread from 0 to 1

#NOISE
plot(grid[,'peeps1'])

#Neighbour list from that plz
grid.sp <- as_Spatial(grid)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(grid.sp, queen=F)

Rcpp::sourceCpp("simfunctions.cpp")

#This is *much* slower! Kind of understandably
x <- optimiseWEIGHTEDAverageAbsoluteContiguousDifference(
  attribute=grid$peeps1, secondpop=grid$peeps2, nblist=neighbours,
  maximise = TRUE, breakval = 150000, cutoff = 0)

grid$opt <- x$attribute
plot(grid[,'opt'])
saveRDS(grid,'local/gridopt1_20x20_maximised.rds')

#Hmm. Not sure what pattern is emerging there. Let's run overnight and see!
#Though computer might turn itself off, but let's see.
x <- tm_shape(grid) +
  tm_fill("opt", n=30, legend.show = F)

tmap_save(tm = x, filename = 'outputs/20x20_DI_uniform_WEIGHTED1_maximised.png', width=10,height=10)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#USE ONLY ONE SIDE PER CELL: OPTIMISE WEIGHTED WITH UNIFORM 50*50 GRID----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ncol = 20
nrow = 20

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")


grid$peeps1 <- runif(nrow(grid))
grid$peeps2 <- 1-grid$peeps1
hist(grid$peeps1)#Evenly spread from 0 to 1

#NOISE
plot(grid[,'peeps1'])

#Neighbour list from that plz
grid.sp <- as_Spatial(grid)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(grid.sp, queen=F)

Rcpp::sourceCpp("simfunctions.cpp")

#This is *much* slower! Kind of understandably
x <- optimiseWEIGHTED_AACD_pick_only_one_sidepercell(
  attribute=grid$peeps1, secondpop=grid$peeps2, nblist=neighbours,
  maximise = F, breakval = 150000, cutoff = 0)

grid$opt <- x$attribute
plot(grid[,'opt'])
saveRDS(grid,'local/gridopt1_20x20_minimised_onlyonesidepercell.rds')

#Hmm. Not sure what pattern is emerging there. Let's run overnight and see!
#Though computer might turn itself off, but let's see.
x <- tm_shape(grid) +
  tm_fill("opt", n=30, legend.show = F)

tmap_save(tm = x, filename = 'outputs/gridopt1_20x20_minimised_onlyonesidepercell.png', width=10,height=10)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CREATE LARGER VERSIONS OF LITTLE ABSTRACT CITIES----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Plan: create original, just scale up
ncol = 7
nrow = 8

#Set all this up with the same geography.
#We have two versions to do.
gridsmooth <- raster(nrow = nrow, ncol = ncol)

gridedges <- raster(nrow = nrow, ncol = ncol) 

values(gridsmooth) <- c(0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.6,0.6,0.05,0.2,0.4,0.4,0.4,0.4,0.4,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
values(gridedges) <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)

#tick
plot(gridedges)
plot(gridsmooth)

#Ah ha
#https://stackoverflow.com/questions/32278825/how-to-change-the-resolution-of-a-raster-layer-in-r
#Though seems to be opposite way round to what that says
gridsmooth <- disaggregate(gridsmooth, fact=7)
gridedges <- disaggregate(gridedges, fact=7)

#So that's not smooth, just larger

gridsmooth <- gridsmooth %>% rasterToPolygons() %>% as("sf")
gridsmooth$id = 1:nrow(gridsmooth)

gridedges <- gridedges %>% rasterToPolygons() %>% as("sf")
gridedges$id = 1:nrow(gridedges)
#tick
plot(gridedges)

gridsmooth <- gridsmooth %>% rename(peeps1 = layer)
gridedges <- gridedges %>% rename(peeps1 = layer)

Rcpp::sourceCpp("simfunctions.cpp")

#Nontorus
smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff=0, torus=F)
edgeAACD <- getAverageAbsoluteContiguousDifference(gridedges$peeps1, ncol = ncol, nrow = nrow,cutoff=0,torus=F)

actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD))

#Randomly permute grid cells, get AACD from permutes.
#Smoothperm. Nice.
#mean is roughly...(same for both cos same cells values just repositioned)
#Torus not huge diff to randomised
mean(getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 2000, cutoff=0, torus=F))

smoothperm = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 10000, cutoff=0, torus=F)
edgeperm = getRepeatedAACDfromPermutedCells(attribute = gridedges$peeps1, ncol = ncol, nrow = nrow, numreps = 10000, cutoff=0, torus=F)

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


#Well that works less well!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#REAL GEOGRAPHY WITH CUTOFF----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Get somewhere. Err. London LSOAs?
lsoas <- st_read('F:/Data/MapPolygons/England/2011/England_lsoa_2011_clipped/england_lsoa_2011_clipped.shp')

#Oh goddam it'll all be London local authority names...
shef <- lsoas %>% 
  filter(grepl(pattern = 'Sheffield',x = .$name))

#Ah, we need an attribute. Did we have some of those elsewhere for Sheffield?
cob <- read_csv('C:/Users/Dan Olner/Dropbox/SheffieldMethodsInstitute/3D_printing2/data/countryOfBirth_Y&H_LSOA.csv')

#345 seems like enough to be going with
#plot(st_geometry(shef))

shef.sp <- as_Spatial(shef)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(shef.sp, queen=F)

neighbours[[25]]

Rcpp::sourceCpp("simfunctions.cpp")
displayAllNeighbours(shef.sp$code,neighbours)


#Via 3D printing script
#We want all of shef if poss, plz
cob <- read_csv('C:/Users/Dan Olner/Dropbox/SheffieldMethodsInstitute/3D_printing2/data/countryOfBirth_Y&H_LSOA.csv')

#Tick
table(shef.sp$code %in% cob$LSOA11CD)

#Subset to Sheffield
cob <- cob[cob$LSOA11CD %in% shef.sp$code,]

#tidy to make easier to read
names(cob) <- gsub('; measures: Value','', names(cob))

#non-white is just 'all categories' minus UK (include all UK)
cob$nonUK <- cob$`Country of Birth: All categories: Country of birth` - 
  (cob$`Country of Birth: Europe: United Kingdom: Total` +
     cob$`Country of Birth: Europe: Great Britain not otherwise specified` +
     cob$`Country of Birth: Europe: United Kingdom not otherwise specified`)

#as % of zone pop
cob$nonUKZoneProp <- (cob$nonUK/cob$`Country of Birth: All categories: Country of birth`)*100
hist(cob$nonUKZoneProp)

#merge with geography - keep only the column we want for now.
#cob_geo <- merge(lsoas[,c('code')], cob[,c('LSOA11CD','nonUKZoneProp')], by.x = 'code', by.y = 'LSOA11CD')

#Add to shef.sp
shef.sp.cob <- merge(shef.sp, cob[,c('LSOA11CD','nonUKZoneProp')], by.x = 'code', by.y = 'LSOA11CD')

#Does that look right? Tick.
shef.sf <- st_as_sf(shef.sp.cob)
plot(shef.sf[,'nonUKZoneProp'])


Rcpp::sourceCpp("simfunctions.cpp")
rez <- getNeighbourIndexAACD(shef.sf$nonUKZoneProp, neighbours)[['AACD']]
getNeighbourIndexAACD(sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp),replace=F), neighbours)[['AACD']]


#OK, so just a little quick permute test to get a mean / null dist for that
nullz <- list()
for(i in 1:5000){
  print(i)
  nullz[[length(nullz)+1]] <- getNeighbourIndexAACD(sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp),replace=F), neighbours)
}

nullz <- unlist(nullz)

ggplot(data.frame(means = nullz), aes(x = means)) +
  geom_density() +
  geom_vline(xintercept = rez, colour = 'red')




Rcpp::sourceCpp("simfunctions.cpp")

#Try optimising
minz <- optimiseGetNeighbourIndexAACD(attribute = shef.sf$nonUKZoneProp, nblist = neighbours,maximise = F,
                                      breakval = 200000,cutoff = 0)
maxz <- optimiseGetNeighbourIndexAACD(attribute = shef.sf$nonUKZoneProp, nblist = neighbours,maximise = T,
                                      breakval = 200000,cutoff = 0)

minz_from_random <- optimiseGetNeighbourIndexAACD(attribute = sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp)), 
                                                  nblist = neighbours,maximise = F,
                                      breakval = 200000,cutoff = 0)

#Maximum reached quite quickly
maxz_from_random <- optimiseGetNeighbourIndexAACD(attribute = sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp)), 
                                                  nblist = neighbours,maximise = T,
                                      breakval = 200000,cutoff = 0)


shef.sf$minz <- minz
shef.sf$maxz <- maxz
shef.sf$minz_from_random <- minz_from_random
shef.sf$maxz_from_random <- maxz_from_random
plot(shef.sf[,c('minz')])
plot(shef.sf[,c('maxz')])
plot(shef.sf[,c('minz_from_random')])
plot(shef.sf[,c('maxz_from_random')])

#The west tends to end up high I suspect cos it has a lot of borders - 
#that maximises ACD if it's v different from surrounding.
getNeighbourIndexAACD(shef.sf$nonUKZoneProp, neighbours)[['AACD']]
#Mean of randomised
getNeighbourIndexAACD(sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp),replace=F), neighbours)[['AACD']]
getNeighbourIndexAACD(shef.sf$minz, neighbours)[['AACD']]
getNeighbourIndexAACD(shef.sf$maxz, neighbours)[['AACD']]




#LOOK AT all pair ACD values to think about cutoff
allz <- getNeighbourIndexAACD(shef.sf$nonUKZoneProp, neighbours)[['allACDs']]

ggplot(data.frame(allz=allz), aes(x=allz))+
  geom_density()

#Ah now that's kind of interesting. Compare to min and max
mindist <- getNeighbourIndexAACD(shef.sf$minz, neighbours)[['allACDs']]
maxdist <- getNeighbourIndexAACD(shef.sf$maxz, neighbours)[['allACDs']]

ggplot(data.frame(ACD=mindist), aes(x=ACD))+
  geom_density()


combo <- bind_rows(
  data.frame(ACD = allz) %>% mutate(type = 'actual'),
  data.frame(ACD = mindist) %>% mutate(type = 'min'),
  data.frame(ACD = maxdist) %>% mutate(type = 'max')
)


ggplot(combo, aes(x=ACD, colour = type))+
  geom_density()

# ggplot(combo, aes(x=ACD, fill = type))+
#   geom_histogram(position="dodge")



#What does dist of ACDs look like for pretend city?
Rcpp::sourceCpp("simfunctions.cpp")
#(Already loaded in above)
##Get neighbour list
gridedges.sp <- as_Spatial(gridedges)
gridsmooth.sp <- as_Spatial(gridsmooth)

#Rook contig plz! We only want bordering cells
gridedges.neighbours <- poly2nb(gridedges.sp, queen=F)
gridsmooth.neighbours <- poly2nb(gridsmooth.sp, queen=F)

#WHUT? Why are they different lengths?? Oh: cutoff should include 0==0. Drop any BELOW.
#length(getNeighbourIndexAACD(attribute = gridsmooth$peeps1,nblist = gridsmooth.neighbours)[['allACDs']])
#length(getNeighbourIndexAACD(attribute = gridedges$peeps1,nblist = gridedges.neighbours)[['allACDs']])
smoothedge <- bind_rows(
  data.frame(ACDs=getNeighbourIndexAACD(attribute = gridsmooth$peeps1,nblist = gridsmooth.neighbours)[['allACDs']],
             type="smooth"),
  data.frame(ACDs=getNeighbourIndexAACD(attribute = gridedges$peeps1,nblist = gridedges.neighbours)[['allACDs']],
             type="edge")
)

ggplot(smoothedge, aes(x=ACDs,fill=type)) +
  geom_histogram(position = "dodge",bins = 8)





#TRY CUTOFF VALUE FOR SHEFFIELD, SEE WHAT THAT DOES TO POSITION OF ACTUAL VS NULL
#Reminder: distribution of actual:
ggplot(data.frame(allz=allz), aes(x=allz))+
  geom_density()


actualresults = list()
permresults = list()

#So if we try cutoffs from 0 to 10?
for(cutoff in seq(from = 0, to = 10, by = 0.25)){
  
  print(cutoff)
  
  nullz <- list()
  for(i in 1:5000){
    nullz[[length(nullz)+1]] <- getNeighbourIndexAACD(sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp),replace=F), neighbours,cutoff = cutoff)[['AACD']]
  }
  
  nullz <- unlist(nullz)
  
  mean_n_quantz <- data.frame(
    mean = mean(nullz), sd = sd(nullz), `2.5%` = quantile(nullz,0.025), `97.5%` = quantile(nullz,0.975),
    cutoff = cutoff
  )
    
  permresults[[length(permresults)+1]] <- mean_n_quantz
  
  #Actual result
  actualresults[[length(actualresults)+1]] <- data.frame(
    AACD = getNeighbourIndexAACD(shef.sf$nonUKZoneProp, neighbours, cutoff = cutoff)[['AACD']], 
    cutoff = cutoff)
  
}


actualallz <- bind_rows(actualresults)
permallz <- bind_rows(permresults)

ggplot(permallz, aes(x = cutoff, y = mean)) +
  geom_point() +
  geom_errorbar(aes(x = cutoff, ymin = `X2.5.`, ymax = `X97.5.`)) +
  xlab('cutoff') +
  ylab('AACD permute: mean & 95% quantiles vs actual AACD') +
  # coord_cartesian(ylim = c(0,1.5)) +
  geom_point(data = actualallz, aes(x = cutoff, y = AACD), size = 2, colour = 'blue')

ggsave('outputs/SHEFFIELD_NONUK_w_actualAACDvalues_cutoffs_notorus.png', width = 12, height = 4)












nullz <- list()
for(i in 1:5000){
  print(i)
  nullz[[length(nullz)+1]] <- getNeighbourIndexAACD(sample(shef.sf$nonUKZoneProp,length(shef.sf$nonUKZoneProp),replace=F), neighbours)
}

nullz <- unlist(nullz)

ggplot(data.frame(means = nullz), aes(x = means)) +
  geom_density() +
  geom_vline(xintercept = rez, colour = 'red')





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SELECTIVE NEIGHBOUR INDEX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)


#Get neighbour list
gridedges.sp <- as_Spatial(gridedges)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(gridedges.sp, queen=F)

#Tick, max four
neighbours[[25]]


#testing... max of two neighbour pairs, that's right
Rcpp::sourceCpp("simfunctions.cpp")
x <- mean(selectiveNeighbourIndex(gridedges$peeps1,neighbours))
y <- mean(selectiveNeighbourIndex(gridsmooth$peeps1,neighbours))
x <- mean(selectiveNeighbourIndex(gridedges$peeps1,neighbours, useNeighMax = F))
y <- mean(selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, useNeighMax = F))


#Compare to null, plot
nullz <- list()
for(i in 1:10000){
  nullz[[length(nullz)+1]] <- mean(selectiveNeighbourIndex(sample(gridedges$peeps1,length(gridedges$peeps1),replace = F),
                                                          neighbours,useNeighMax = F))
}

nullz <- unlist(nullz)

ggplot(data.frame(null = nullz), aes(x=nullz)) +
  geom_density() +
  geom_vline(xintercept = mean(x), colour='red') +
  geom_vline(xintercept = mean(y), colour='green')


#Well that's no different, huh?
ggsave('outputs/selectiveneighbourindexnull.png',width=7,height=7)


#That pattern looks... blocky! We shall see. Now for optimising a random space
ncol = 20
nrow = 20

#Set all this up with the same geography.
grid <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")


grid$peeps1 <- runif(nrow(grid))
grid$peeps2 <- 1-grid$peeps1
hist(grid$peeps1)#Evenly spread from 0 to 1

#NOISE
plot(grid[,'peeps1'])

#Neighbour list from that plz
grid.sp <- as_Spatial(grid)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(grid.sp, queen=F)




minz <- optimiseSelectiveNeighbourIndex(
  attribute = grid$peeps1, nblist = neighbours,
  maximise = F, breakval = 50000,useNeighMax = T)

grid$minz <- minz
plot(grid[,'minz'])

grid$maxz <- optimiseSelectiveNeighbourIndex(
  attribute = grid$peeps1, nblist = neighbours,
  maximise = T, breakval = 50000,useNeighMax = T)
plot(grid[,'maxz'])




#TEST SAME WITH THRESHOLD (default is useNeighMax = T)
Rcpp::sourceCpp("simfunctions.cpp")

#Check range to see what good threshold might be 
range(selectiveNeighbourIndex(gridedges$peeps1,neighbours, threshold = 0))
range(selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, threshold = 0))
selectiveNeighbourIndex(gridedges$peeps1,neighbours, threshold = 0)
selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, threshold = 0)

x <- mean(selectiveNeighbourIndex(gridedges$peeps1,neighbours, threshold = 0))
y <- mean(selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, threshold = 0))

#Compare to null, plot
nullz <- list()
for(i in 1:10000){
  nullz[[length(nullz)+1]] <- mean(selectiveNeighbourIndex(sample(gridedges$peeps1,length(gridedges$peeps1),replace = F),
                                                           neighbours,useNeighMax = F))
}

nullz <- unlist(nullz)

ggplot(data.frame(null = nullz), aes(x=nullz)) +
  geom_density() +
  geom_vline(xintercept = mean(x), colour='red') +
  geom_vline(xintercept = mean(y), colour='green')




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#BRINGING WEIGHTED NEIGHBOUR INDEX BACK INTO R----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#From C++. For ease of working with / adding new shizzle.
#I already have code for it, used for debugging. But it doesn't yet (I don't think)
#Find any actual values. Let's looksee.

#See notes section "new measure?"
#Use G's little sample cities again - will maybe scale them up at some point
# ncol = 7
# nrow = 8
# 
# #Set all this up with the same geography.
# #We have two versions to do.
# gridsmooth <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")
# gridsmooth$id = 1:nrow(gridsmooth)
# 
# gridedges <- raster(nrow = nrow, ncol = ncol) %>% rasterToPolygons() %>% as("sf")
# gridedges$id = 1:nrow(gridedges)


#RASTERTOPOLYGONS seems to return a rectangular grid. 
#We want exactly square to keep distance calcs correct.
#Use one I made earlier...


gridsmooth <- st_read('saves/shapefiles/7by8gridSquare.shp')
gridedges <- st_read('saves/shapefiles/7by8gridSquare.shp')

#For other work
grid.rci <- st_read('saves/shapefiles/7by8gridSquare.shp')

#Of course, it would be oriented completely differently. Can I use ID to fix?
plot(gridsmooth)

#Reorient via ID so order is top-left across
neworder <- list()
for(i in 1:7) neworder[[length(neworder)+1]] <- seq(from=i, to=56, by=7)

gridsmooth$id <- unlist(neworder)
gridsmooth <- gridsmooth %>% arrange(id)

gridedges$id <- unlist(neworder)
gridedges <- gridedges %>% arrange(id)

grid.rci$id <- unlist(neworder)
grid.rci <- grid.rci %>% arrange(id)

#Keep only 7*7 grid
grid.rci <- grid.rci[1:49,]
st_write(grid.rci,'grid.rci.shp',delete_layer = T)


  
#Coding values directly from G (left to right)
#G uses single value - for DI for two pops, the other peeps are 1-x
gridsmooth$peeps1 <- c(0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.8,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.8,0.8,0.05,0.2,0.4,0.6,0.6,0.6,0.6,0.05,0.2,0.4,0.4,0.4,0.4,0.4,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.2,0.2,0.2,0.2,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
gridsmooth$peeps2 <- 1 - gridsmooth$peeps1

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,0.2,0.2,0.2,0.05,0.2,0.2,0.4,0.2,0.2,0.05,0.8,0.05,0.6,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.8,0.8,0.8,0.05,0.4,0.2,0.05,0.05,0.8,0.05,0.05,0.6,0.2,0.6,0.05,0.05,0.05,0.6,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1- gridedges$peeps1

plot(gridsmooth)
plot(gridedges)

#Check order is: start top-left, go left to right. Tick.
# gridsmooth$peeps1[8] <- 100
# plot(gridsmooth)

#Get neighbour list
gridedges.sp <- as_Spatial(gridedges)

#Rook contig plz! We only want bordering cells
neighbours <- poly2nb(gridedges.sp, queen=F)

#Tick, max four
neighbours[[25]]

displayAllNeighbours(gridedges$peeps1,neighbours)



#SHIFTED NEMA WEIGHTED BORDER CALC TO FUNCTION
weightedAACD(gridedges$peeps1,neighbours)


#Null dist
y = list()

for(i in 1:1000){
  
  cat(i,'\n')
  
  y[[length(y)+1]] <- weightedAACD(sample(gridedges$peeps1,replace = F),neighbours)
  # y <- weightedAACD(sample(gridedges$peeps1,replace = F), neighbours)
  
}

plot(density(unlist(y)))


#Routine for optimising. Permute values, seek max/min
x = proc.time()

var <- gridedges$peeps1

for(i in 1:10000){
  
  y <- weightedAACD(var, neighbours)
  
  newvar = var
  
  #Swap single cell
  # swap1 = var[ as.integer(runif(1)*length(var))  ]
  swap1 = sample(1:length(newvar),1)
  swap2 = sample(1:length(newvar),1)
  # swap1 = newvar[ sample(1:length(newvar),1)  ]
  # swap2 = newvar[ sample(1:length(newvar),1)  ]
  
  #Swap em
  temp = newvar[swap1]
  newvar[swap1] = newvar[swap2]
  newvar[swap2] = temp
  
  z = weightedAACD(newvar, neighbours)
  
  #Minimise
  if(z < y) var = newvar
  #Maximise
  # if(z > y) var = newvar
  
}

proc.time()-x

gridedges$opt_peeps <- var
plot(gridedges)

#Speed comparison to C++ plz?
x = proc.time()

minz <- optimiseSelectiveNeighbourIndex(
  attribute = gridedges$peeps1, nblist = neighbours,
  maximise = F, breakval = 10000,useNeighMax = T)

proc.time()-x

#C++ 1.74 seconds
#R 84.25 second
#C++ is 50 times faster
#(For 8*7 grid / 10000 iterations)

#Is C++ and R giving same values? Tick. That seems to be working then...
#[[1]] is index, [[2]] contains actual border values in list
x <- weightedAACD(gridedges$peeps1,neighbours)
x[[1]]
mean(selectiveNeighbourIndex(gridedges$peeps1,neighbours))
x <- weightedAACD(gridedges$opt_peeps,neighbours)
x[[1]]
mean(selectiveNeighbourIndex(gridedges$opt_peeps,neighbours))


#Same actual order of values returned? TICK.
x <- weightedAACD(gridedges$peeps1,neighbours)[[2]]
y <- selectiveNeighbourIndex(gridedges$peeps1,neighbours)

#So actually, what that means:
#I can get vector from C++ and then get another vector for Moran's I
#And combine them.

#~~~~~~~~~~~~~~~~
#LOCAL MORANS----
#~~~~~~~~~~~~~~~~

#Moran's from spdep
x = nb2listw(neighbours)
moran.test(gridedges$peeps1,x)


#So. The plan:
#Find Moran's I on each side of pairs. Add them. Keep as edge index.
#Maintain order run as before

#Ah, actually: to get the full 9 zone neighbour list, easiest to use
#Queen contig
neighbours.queen <- poly2nb(gridedges.sp, queen=T)
gridsmooth.sp <- as_Spatial(gridsmooth)


debugonce(localMoransI_onEachSideOfBorder)

localMoransI_onEachSideOfBorder(gridedges$peeps1,neighbours,neighbours.queen,gridedges.sp)
localMoransI_onEachSideOfBorder(gridsmooth$peeps1,neighbours,neighbours.queen,gridsmooth.sp)


#Test optimising just on Moran's. What does that look like?
#(Very small grid, might be worth trying with larger...)
x = proc.time()

var <- gridedges$peeps1
origspatial <- gridedges.sp#Use this to get local neighbour values for Moran's
#Use copy so we can update this too. Might be an idea just to use that to save confusion, really...

#SLOOOW!
for(i in 1:100){
  
  cat(i,"\n")
  
  y <- localMoransI_onEachSideOfBorder(var, neighbours, neighbours.queen, origspatial)
  
  newvar = var
  newspatial = origspatial
  
  #Swap single cell
  swap1 = sample(1:length(newvar),1)
  swap2 = sample(1:length(newvar),1)
  
  #Swap em (both in var and spatial version... again, should make one copy at some point to keep simpler)
  temp = newvar[swap1]
  newvar[swap1] = newvar[swap2]
  newspatial$peeps1[swap1] = newvar[swap2]
  
  newvar[swap2] = temp
  newspatial$peeps1[swap2] = temp
  
  z <- localMoransI_onEachSideOfBorder(newvar, neighbours, neighbours.queen, newspatial)
  
  #Minimise
  if(z[[1]] > y[[1]]) {
    
    var = newvar
    origspatial = newspatial
    
  }
  #Maximise
  # if(z > y) var = newvar
  
}

proc.time()-x

gridedges$opt_peeps <- var
plot(gridedges)


# newvar = newvar[sample(1:length(newvar), replace = F)]
# #Still same values
# newvar[order(newvar)]==gridedges$peeps1[order(gridedges$peeps1)]
# newvar[order(newvar)]==gridsmooth$peeps1[order(gridsmooth$peeps1)]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TESTING BASIC DECENTRALISATION INDEX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Looking at Massey/denton currently. Let's find a geography to play with. Err...
shp <- st_read('C:/Users/Dan Olner/Dropbox/SheffieldMethodsInstitute/CountryOfBirthOpenDataSets/data/gb_shapefile/gb_altered_wards_n_postcodesectors_w_lookup.shp')

lon.data <- read_csv('C:/Users/Dan Olner/Dropbox/SheffieldMethodsInstitute/CountryOfBirthOpenDataSets/data/countryofbirth/countryOfBirth_GreatBritain_5_Census.csv', guess_max = 15000)

#Nabbed from repo
#library(conflicted)
#conflict_prefer("select", "dplyr") 

props <- lon.data %>%
  mutate(rowsums = rowSums(dplyr::select(.,England:`Rest of world`))) %>% 
  mutate_at(vars(England:`Rest of world`), funs(  ((.)/rowsums)*100  )) %>% 
  dplyr::select(-rowsums)

#Add a column for total of non-UK-born
props <- props %>%
  mutate(nonUK = `Irish Republic`+India+Pakistan+Europe+`Rest of world`,
         UK = England+Scotland+Wales+`Rest of UK`,
        Total = `Irish Republic`+India+Pakistan+Europe+`Rest of world`+England+Scotland+Wales+`Rest of UK`
         )

#Merge spatial zones with CoB data
join <- left_join(shp,props,by='zone')

#Subset to London TTWA
london <- join %>% filter(ttwa == 'London')
plot(st_geometry(london))


#..

#So in theory, didn't need to do proportions anyway, did we? In fact, that was a bad idea...
#Try again!
counts <- lon.data %>%
  mutate(
    count_UK = rowSums(dplyr::select(.,England:Wales,`Rest of UK`)),
    count_nonUK = rowSums(dplyr::select(.,`Irish Republic`,India:`Rest of world`))
         )

join <- left_join(shp,counts,by='zone')

#Subset to London TTWA
london <- join %>% filter(ttwa == 'London')

london11 <- london %>% filter(censusYear==2011) %>% dplyr::select(zone,count_UK,count_nonUK)


#It perhaps has to actually be proportions across zones. Let's see what happens.
#First-up, need to work out order of zones by distance from centre.
#Pick a centre zone... 
#01ALFJ will do.
#is this zero cos they're in contact? Better to use centroids?
# st_distance(london11[1,],london11[2,])
# 
# #Yup.
# london11.centroids <- st_centroid(london11)
# st_distance(london11.centroids[1,],london11.centroids[2,])

#So.
london11$distfromcentre_01ALFJ <- st_distance(london11.centroids[london11.centroids$zone=='01ALFJ',],london11.centroids,by_element = T)

#Sort by distance from centre
london11 <- london11 %>% arrange(distfromcentre_01ALFJ)
#Get rid of geography (or it finds value for the geometry...)
london11.nogeog <- london11 %>% st_set_geometry(NULL)

#RCI
tot = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(london)){
  tot[[length(tot)+1]] <- cumsum(london11.nogeog[1:(i-1),'count_UK'])
    
    
    london11.nogeog[(i-1),'count_UK'] * london11.nogeog[i,'count_nonUK']
}




for(i in 2:nrow(london)){
  tot[[length(tot)+1]] <- london11.nogeog[(i-1),'count_UK'] * london11.nogeog[i,'count_nonUK']
}

tot2 = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(london)){
  tot2[[length(tot)+1]] <- london11.nogeog[(i),'count_UK'] * london11.nogeog[i-1,'count_nonUK']
}

#Sum and substract the two...



#DOING IT WITH LONDON, MAYBE A BAD IDEA. THAT'S A LOTTA ZONES.
counts <- lon.data %>%
  mutate(
    count_UK = rowSums(dplyr::select(.,England:Wales,`Rest of UK`)),
    count_nonUK = rowSums(dplyr::select(.,`Irish Republic`,India:`Rest of world`))
  )

join <- left_join(shp,counts,by='zone')

#Subset to London TTWA
shef <- join %>% filter(grepl(pattern = 'Sheffield', .$ttwa))
shef11 <- shef %>% filter(censusYear==2011) %>% dplyr::select(zone,count_UK,count_nonUK)


#05CGFF is centralish
shef11.centroids <- st_centroid(shef11)
shef11$distfromcentre_01ALFJ <- st_distance(shef11.centroids[shef11.centroids$zone=='05CGFF',],shef11.centroids,by_element = T)

#Sort by distance from centre
shef11 <- shef11 %>% arrange(distfromcentre_01ALFJ)
#Get rid of geography (or it finds value for the geometry...)
shef11.nogeog <- shef11 %>% st_set_geometry(NULL)

#cumsum(shef11.nogeog$count_UK)

#RCI
tot = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(shef11.nogeog)){
  
  x = cumsum(shef11.nogeog[1:(i-1),'count_UK'])
  x = x[length(x)]#get final cumulative sum
  
  y = cumsum(shef11.nogeog[2:(i),'count_nonUK'])
  y = y[length(y)]#get final cumulative sum
  
    tot[[length(tot)+1]] <- x*y
  
}

#Sum all those results
tot = sum(unlist(tot))

tot2 = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(shef11.nogeog)){
  
  x = cumsum(shef11.nogeog[2:(i),'count_UK'])
  x = x[length(x)]#get final cumulative sum
  
  y = cumsum(shef11.nogeog[1:(i-1),'count_nonUK'])
  y = y[length(y)]#get final cumulative sum
  
  tot2[[length(tot2)+1]] <- x*y
  
}

#Sum all those results
tot2 = sum(unlist(tot2))



#OK, THINK EACH NEEDS TO BE PROPORTION FOR EACH COLUMN
shef11.nogeog <- shef11.nogeog %>% 
  mutate_at(vars(count_UK:count_nonUK), funs(prop = ./sum(.)))

#sum(shef11.nogeog$count_nonUK_prop)#tick


tot = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(shef11.nogeog)){
  
  x = cumsum(shef11.nogeog[1:(i-1),'count_UK_prop'])
  x = x[length(x)]#get final cumulative sum
  
  y = cumsum(shef11.nogeog[1:(i),'count_nonUK_prop'])
  y = y[length(y)]#get final cumulative sum
  
    tot[[length(tot)+1]] <- x*y
  
}

#Sum all those results
tot = sum(unlist(tot))

tot2 = list()

#Sum cum proportions multiplied for first lot
for(i in 2:nrow(shef11.nogeog)){
  
  x = cumsum(shef11.nogeog[1:(i),'count_UK_prop'])
  x = x[length(x)]#get final cumulative sum
  
  y = cumsum(shef11.nogeog[1:(i-1),'count_nonUK_prop'])
  y = y[length(y)]#get final cumulative sum
  
  tot2[[length(tot2)+1]] <- x*y
  
}

#Sum all those results
tot2 = sum(unlist(tot2))

tot-tot2

#Negative: y members are closer to centre, relatively... Tick.
plot(shef11[,2:3])



#function version
#ASSUMES ALREADY ORDERED BY DISTANCE
relativeCentralisationIndex(shef11.nogeog,'count_UK_prop','count_nonUK_prop')
relativeCentralisationIndex(shef11.nogeog,'count_nonUK_prop','count_UK_prop')


#Check against Meng Le's version
rci(x = shef11.nogeog$count_UK,y = shef11.nogeog$count_nonUK,sort.var = shef11.nogeog$distfromcentre_01ALFJ)
rci(x = shef11.nogeog$count_nonUK,y = shef11.nogeog$count_UK,sort.var = shef11.nogeog$distfromcentre_01ALFJ)

#check speed diff
x <- proc.time()
replicate(1000,relativeCentralisationIndex(shef11.nogeog,'count_UK_prop','count_nonUK_prop'))
proc.time()-x

#25 times faster!
x <- proc.time()
replicate(1000,rci(x = shef11.nogeog$count_UK,y = shef11.nogeog$count_nonUK,sort.var = shef11.nogeog$distfromcentre_01ALFJ))
proc.time()-x

#"Positive values indicate that X members are located closer to the city centre"... Massey / Denton
#"Positive values imply that poor individuals are relatively more concentrated near the centre compared with non-poor" Meng Le / Gwilym (where "poor" are equiv to x)

#Let's change names to x and y to make clear
df <- shef11.nogeog %>% 
  dplyr::select(x = count_UK, y = count_nonUK, distance = distfromcentre_01ALFJ)

df.sp <- shef11 %>% 
  dplyr::select(x = count_UK, y = count_nonUK, distance = distfromcentre_01ALFJ)


rci(x = df$x,y = df$y,sort.var = df$distance)

plot(df.sp[,c('x','y')])

st_write(df.sp,'saves/shapefiles/check.shp', delete_layer = T)


#Or using the GP/MLZ paper notation:
df <- shef11.nogeog %>% 
  dplyr::select(a = count_UK, b = count_nonUK, distance = distfromcentre_01ALFJ, zone)

df.sp <- shef11 %>% 
  dplyr::select(a = count_UK, b = count_nonUK, distance = distfromcentre_01ALFJ, zone)


rci_ab(a = df$a,b = df$b,sort.var = df$distance)

plot(df.sp[,c('a','b')])


#Edited version to make "a" definitely in the outside ring
df$a <- 2

#Got outer edge from QGIS
df$a[df$zone %in% c(
'05CGGD',
'05CGGE',
'05CGFG',
'05CFFE',
'05CFFR',
'05CFFJ',
'05CFFD',
'05CFFP',
'05CFFT',
'05CFFA',
'05CFFN',
'18FTFL',
'18FTFK',
'18FTFX',
'18FTFF',
'18FTFG',
'18FTFB',
'18FQFN',
'18FQFQ',
'05CGFM'
)] <- 10

df.sp$a <- 2

#Got outer edge from QGIS
df.sp$a[df.sp$zone %in% c(
'05CGGD',
'05CGGE',
'05CGFG',
'05CFFE',
'05CFFR',
'05CFFJ',
'05CFFD',
'05CFFP',
'05CFFT',
'05CFFA',
'05CFFN',
'18FTFL',
'18FTFK',
'18FTFX',
'18FTFF',
'18FTFG',
'18FTFB',
'18FQFN',
'18FQFQ',
'05CGFM'
)] <- 10


plot(df.sp[,c('a','b')])

rci_ab(a = df$a,b = df$b,sort.var = df$distance)

st_write(df.sp,'saves/shapefiles/check_altered.shp', delete_layer = T)


#~~~~~~~~~~~~~~~~~~~~~~
#TESTING: ORDERING GRID ZONES BY DISTANCE; GROUPING ONES THAT ARE THE SAME----
#~~~~~~~~~~~~~~~~~~~~~~

#ie. rook contig cells will be grouped
#And queen contig a different group
#If we go out by more than one order, centroid distance then grouping should still work...

#With new shapefile with actual squares, not rectangles

#Apparently... Oh except, can't do on just one dimension?
#https://r-spatial.github.io/sf/articles/sf3.html
# gridedges.transform <- gridedges * 0.5


#Test example (using prev. run grid creation code above.)
#Get neighbour list
gridedges.sp <- as_Spatial(gridedges)

#Rook contig plz! We only want bordering cells
neighbours.queen <- poly2nb(gridedges.sp, queen=T)



#Pick an example with eight neighbours
neighbours.queen[[25]]

#Get those simple-features cells: 9 of em
littlegroup <- gridedges[c(25,neighbours.queen[[25]]),]

#Centroids
littlegroup.c <- st_centroid(littlegroup)

#Distance from central point. First row is centre, selected that above.
littlegroup$distances <- st_distance(littlegroup.c[1,],littlegroup.c[1:9,], by_element = T)

#Those should be in three distance groups: centre, rook, queen... now we have squares, tick!
unique(littlegroup$distances)


#So. If we want three groups for our RCI,
#Need values to be ... well, population per unit of area probably, right?
#I think sf might take care of that for us...
littlegroup.rci <- littlegroup %>%
  group_by(distances) %>% 
  summarise(peeps1 = mean(peeps1), distance = mean(distances))

plot(littlegroup.rci)

#Is mean of four rook zones 0.8...? Yup
tmap_mode("view")
qtm(littlegroup %>% dplyr::select(peeps1))
qtm(littlegroup.rci %>% dplyr::select(peeps1,distance))

#Y is "first"
rci(x = 1-littlegroup.rci$peeps1, y = littlegroup.rci$peeps1,sort.var = littlegroup.rci$distance)



#Just also checking: for gridedges as a whole, peeps1 is more centralised, right?
gridedges.centr <- st_centroid(gridedges)

#Distance from central point. First row is centre, selected that above.
gridedges$distances <- st_distance(gridedges.centr[25,],gridedges.centr, by_element = T)

#Yup, positive, just.
rci(1-gridedges$peeps1,gridedges$peeps1,gridedges$distances)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TEST RUNNING RCI ON EACH SIDE OF BORDER, FINDING VALUE FOR PEACHY INDEX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#What shall we pass in? Something easier to work with
#Just pass in the full sf?
# gridedges.sp <- as_Spatial(gridedges)
# 
# #Rook contig plz! We only want bordering cells
# neighbours <- poly2nb(gridedges.sp, queen=F)
# #Except we also want to get little group surrounding cell...
# neighbours.queen <- poly2nb(gridedges.sp, queen=T)

debugonce(centralisationIndex_onEachSideOfBorder)

x <- centralisationIndex_onEachSideOfBorder(data.sf = gridedges, var = "peeps1")

#Get absolute contig diff values (via Nema code in permutation_tests)
#Shifted over to function
debugonce(aacd)
aacd(data.sf = gridedges,var = "peeps1")


#~~~~~~~~~~~~~~~~~~~~~~~~~~
#Exponenting? (Answer: newp!)----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

#Just testing whether exponenting values (so larger become much larger) helps in any way?
#Note: exponentiating below zero values will shrink em. Not sure that gets us anywhere...
#Oh no, it still does, curve-wise... (Lowers low values rather than raising high values, same curve)
x = seq(from = 0.01, to = 1, by = 0.01)
y = data.frame(id = x, x = x, y = x^2, z = x^3, a = x^4)

ggplot(y %>% gather(key=key,value=value,x:a), aes(x = id, y = value, colour = key)) +
  geom_line()


y$diff = y$z-y$x
plot(y$diff)

#Trying...
expvalue = 5

gridedges.exp <- gridedges %>% 
  mutate(peeps1exp = peeps1^expvalue)
gridsmooth.exp <- gridsmooth %>% 
  mutate(peeps1exp = peeps1^expvalue)

x <- aacd(data.sf = gridedges.exp,var = "peeps1exp")[[1]]
y <- aacd(data.sf = gridsmooth.exp,var = "peeps1exp")[[1]]

#Get null. Permute exponentiated gridedges peeps1 values
nullz <- list()
for(i in 1:1000){
  cat(i,'\n')
  
  z <- aacd(
    gridedges.exp %>% mutate(permute = sample(peeps1exp,length(peeps1exp),replace = F)),
    var = "permute")
  
  nullz[[length(nullz)+1]] <- z$mean
  
}

nullz <- unlist(nullz)

ggplot(data.frame(null = nullz), aes(x=nullz)) +
  geom_density() +
  geom_vline(xintercept = x, colour='red') +
  geom_vline(xintercept = y, colour='green')


#OR INSTEAD TRY EXPONENTING IN THE ACTUAL FUNCTION ON EACH AACD RESULT
#Not sure that can make much difference, but... NEWP! Thought it did for a second there.
expz = 1

# debugonce(aacd)
x <- aacd(data.sf = gridedge_centre,var = "peeps1",expz = expz)[[1]]
y <- aacd(data.sf = gridsmooth,var = "peeps1",expz = expz)[[1]]

nullz <- list()
for(i in 1:50){
  cat(i,'\n')
  
  z <- aacd(
    gridedge_centre %>% mutate(permute = sample(peeps1,length(peeps1),replace = F)),
    var = "permute", expz = expz)
  
  nullz[[length(nullz)+1]] <- z$mean
}

nullz <- unlist(nullz)

ggplot(data.frame(null = nullz), aes(x=nullz)) +
  geom_density() +
  geom_vline(xintercept = x, colour='red') +
  geom_vline(xintercept = y, colour='green')



#~~~~~~~~~~~~~~~~~~~~~~~~~~
#RCI * AACD combo----
#~~~~~~~~~~~~~~~~~~~~~~~~~~

#Check on gridedges, gridsmooth

#Get null first. Slow, let's stick to 100 for quick test.
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

x <- proc.time()
foreach(icount(100),.combine = c) %dopar% {
  
  cat('ping!\n')
  
  library(sf)
  library(spdep)
  
  permuted <- gridedges
  permuted$permute <- sample(permuted$peeps1,length(permuted$peeps1),replace = F)
  
  a <- aacd(permuted, var = "permute")
  
  b <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "permute")
  
  #Multiply those two
  mean(a$values * b$values)

}
proc.time()-x



x <- proc.time()
rci.aacd.nullz <- list()
for(i in 1:100){
  # cat(i,'\n')
  
  #Reminder: edges, smooth - same cells, just in different positions
  # permuted <- gridedges %>% mutate(permute = sample(peeps1,length(peeps1),replace = F))
  permuted <- gridedges
  permuted$permute <- sample(permuted$peeps1,length(permuted$peeps1),replace = F)
  
  
  a <- aacd(permuted, var = "permute")
  
  b <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "permute")
  
  #Multiply those two
  rci.aacd.nullz[[length(rci.aacd.nullz)+1]] <- mean(a$values * b$values)
  
}
proc.time()-x

#Parallel is three times faster (though we can't see output...)

rci.aacd.nullz <- unlist(rci.aacd.nullz)

plot(density(rci.aacd.nullz))


#Get same result for edges and smooth
a.edges <- aacd(gridedges, var = "peeps1")
b.edges <- centralisationIndex_onEachSideOfBorder(gridedges, var = "peeps1")

#Multiply those two
edges.rci.aacd <- mean(a.edges$values * b.edges$values)


a.smooth <- aacd(gridsmooth, var = "peeps1")
b.smooth <- centralisationIndex_onEachSideOfBorder(gridsmooth, var = "peeps1")

#Multiply those two
smooth.rci.aacd <- mean(a.smooth$values * b.smooth$values)


ggplot(data.frame(null = rci.aacd.nullz), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.rci.aacd, colour='red') +
  geom_vline(xintercept = smooth.rci.aacd, colour='green')



#~~~~~~~~~~~~~~~~~~~~
#RCI * selective neighbour index (max) combo----
#~~~~~~~~~~~~~~~~~~~~

#Get null
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

x <- proc.time()

#Doesn't like C++
# null.rci.neigh <- foreach(icount(10),.combine = c) %dopar% {
null.rci.neigh <- sapply(1:100, function(x) {

  cat(x,'\n')
  
  # library(sf)
  # library(spdep)
  # Rcpp::sourceCpp("simfunctions.cpp")
  
  permuted <- gridedges
  permuted$permute <- sample(permuted$peeps1,length(permuted$peeps1),replace = F)
  
  a <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "permute")
  
  b <- selectiveNeighbourIndex(permuted$permute,neighbours)#Neighbours should still be correct...
  
  #Multiply those two
  mean(a$values * b)
  
})

# }

proc.time()-x

plot(density(null.rci.neigh))


a.edges <- centralisationIndex_onEachSideOfBorder(gridedges, var = "peeps1")
b.edges <- selectiveNeighbourIndex(gridedges$peeps1, neighbours)

#Multiply those two
edges.rci.neigh <- mean(a.edges$values * b.edges)


a.smooth <- centralisationIndex_onEachSideOfBorder(gridsmooth, var = "peeps1")
b.smooth <- selectiveNeighbourIndex(gridsmooth$peeps1, neighbours)

#Multiply those two
smooth.rci.neigh <- mean(a.smooth$values * b.smooth)


ggplot(data.frame(null = null.rci.neigh), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.rci.neigh, colour='red') +
  geom_vline(xintercept = smooth.rci.neigh, colour='green')




#~~~~~~~~~~~~~~~~~~~~
#RCI * selective neighbour index (min) combo----
#~~~~~~~~~~~~~~~~~~~~

#Get null
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

x <- proc.time()

#Doesn't like C++
# null.rci.neigh <- foreach(icount(10),.combine = c) %dopar% {
null.rci.neigh.min <- sapply(1:100, function(x) {

  cat(x,'\n')
  
  # library(sf)
  # library(spdep)
  # Rcpp::sourceCpp("simfunctions.cpp")
  
  permuted <- gridedges
  permuted$permute <- sample(permuted$peeps1,length(permuted$peeps1),replace = F)
  
  a <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "permute")
  
  b <- selectiveNeighbourIndex(permuted$permute,neighbours, useNeighMax = F)#Neighbours should still be correct...
  
  #Multiply those two
  mean(a$values * b)
  
})

# }

proc.time()-x

plot(density(null.rci.neigh))


a.edges <- centralisationIndex_onEachSideOfBorder(gridedges, var = "peeps1")
b.edges <- selectiveNeighbourIndex(gridedges$peeps1, neighbours, useNeighMax = F)

#Multiply those two
edges.rci.neigh.min <- mean(a.edges$values * b.edges)


a.smooth <- centralisationIndex_onEachSideOfBorder(gridsmooth, var = "peeps1")
b.smooth <- selectiveNeighbourIndex(gridsmooth$peeps1, neighbours, useNeighMax = F)

#Multiply those two
smooth.rci.neigh.min <- mean(a.smooth$values * b.smooth)


ggplot(data.frame(null = null.rci.neigh.min), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.rci.neigh.min, colour='red') +
  geom_vline(xintercept = smooth.rci.neigh.min, colour='green')




#~~~~~~~~~~~~~~~~~~~~
#MORAN * AACD----
#~~~~~~~~~~~~~~~~~~~~

#Get null
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

x <- proc.time()

#Doesn't like C++
# null.rci.neigh <- foreach(icount(10),.combine = c) %dopar% {

null.moran.aacd <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  # library(sf)
  # library(spdep)
  # Rcpp::sourceCpp("simfunctions.cpp")
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  permuted.sp <- as_Spatial(permuted)
  
  
  a <- aacd(permuted, var = "peeps1")
  
  #Urgh, this has hardcoded variable names! Needs to be peeps1. What on earth was I thinking?
  # debugonce(localMoransI_onEachSideOfBorder)
  b <- localMoransI_onEachSideOfBorder(permuted$peeps1, neighbours, neighbours.queen, permuted.sp)
  
  #Multiply those two
  mean(a$values * b$values)
  
})

# }

proc.time()-x

plot(density(null.moran.aacd))


a.edges <- aacd(gridedges, var = "peeps1")
b.edges <- localMoransI_onEachSideOfBorder(gridedges$peeps1, neighbours, neighbours.queen, gridedges.sp)

#Multiply those two
edges.moran.aacd <- mean(a.edges$values * b.edges$values)


a.smooth <- aacd(gridsmooth, var = "peeps1")
b.smooth <- localMoransI_onEachSideOfBorder(gridsmooth$peeps1, neighbours, neighbours.queen, gridsmooth.sp)

#Multiply those two
smooth.moran.aacd <- mean(a.smooth$values * b.smooth$values)


ggplot(data.frame(null = null.moran.aacd), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.moran.aacd, colour='red') +
  geom_vline(xintercept = smooth.moran.aacd, colour='green')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#JUST PLAIN MORAN BY ITSELF....----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Nowt unusual there

x <- proc.time()

null.moran <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  permuted.sp <- as_Spatial(permuted)
  
  
  b <- localMoransI_onEachSideOfBorder(permuted$peeps1, neighbours, neighbours.queen, permuted.sp)
  
  #Multiply those two
  mean(b$values)
  
})

# }

proc.time()-x

plot(density(null.moran))


b.edges <- localMoransI_onEachSideOfBorder(gridedges$peeps1, neighbours, neighbours.queen, gridedges.sp)

#Multiply those two
edges.moran <- mean(b.edges$values)


b.smooth <- localMoransI_onEachSideOfBorder(gridsmooth$peeps1, neighbours, neighbours.queen, gridsmooth.sp)

#Multiply those two
smooth.moran <- mean(b.smooth$values)


ggplot(data.frame(null = null.moran), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.moran, colour='red') +
  geom_vline(xintercept = smooth.moran, colour='green')




#~~~~~~~~~~~~~~~~~~~~
#MORAN * SELECTIVE NEIGHBOUR (MAX)----
#~~~~~~~~~~~~~~~~~~~~

#Get null
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

x <- proc.time()

#Doesn't like C++
# null.rci.neigh <- foreach(icount(10),.combine = c) %dopar% {

null.moran.neigh <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  # library(sf)
  # library(spdep)
  # Rcpp::sourceCpp("simfunctions.cpp")
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  permuted.sp <- as_Spatial(permuted)
  
  
  a <- selectiveNeighbourIndex(permuted$peeps1,neighbours, useNeighMax = T)#Neighbours should still be correct...
  
  #Urgh, this has hardcoded variable names! Needs to be peeps1. What on earth was I thinking?
  # debugonce(localMoransI_onEachSideOfBorder)
  b <- localMoransI_onEachSideOfBorder(permuted$peeps1, neighbours, neighbours.queen, permuted.sp)
  
  #Multiply those two
  mean(a * b$values)
  
})

# }

proc.time()-x

plot(density(null.moran.neigh))


a.edges <- selectiveNeighbourIndex(gridedges$peeps1, neighbours)
b.edges <- localMoransI_onEachSideOfBorder(gridedges$peeps1, neighbours, neighbours.queen, gridedges.sp)

#Multiply those two
edges.moran.neigh <- mean(a.edges * b.edges$values)


a.smooth <- selectiveNeighbourIndex(gridsmooth$peeps1, neighbours)
b.smooth <- localMoransI_onEachSideOfBorder(gridsmooth$peeps1, neighbours, neighbours.queen, gridsmooth.sp)

#Multiply those two
smooth.moran.neigh <- mean(a.smooth * b.smooth$values)


ggplot(data.frame(null = null.moran.neigh), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.moran.neigh, colour='red') +
  geom_vline(xintercept = smooth.moran.neigh, colour='green')






#~~~~~~~~~~~~~~~~~~~~
#MORAN * SELECTIVE NEIGHBOUR (MIN)----
#~~~~~~~~~~~~~~~~~~~~

#Get null
# library(doParallel)
# cl <- makeCluster(4)
# registerDoParallel(cl)

x <- proc.time()

#Doesn't like C++
# null.rci.neigh <- foreach(icount(10),.combine = c) %dopar% {

null.moran.neigh.min <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  # library(sf)
  # library(spdep)
  # Rcpp::sourceCpp("simfunctions.cpp")
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  permuted.sp <- as_Spatial(permuted)
  
  
  a <- selectiveNeighbourIndex(permuted$peeps1,neighbours, useNeighMax = F)#Neighbours should still be correct...
  
  #Urgh, this has hardcoded variable names! Needs to be peeps1. What on earth was I thinking?
  # debugonce(localMoransI_onEachSideOfBorder)
  b <- localMoransI_onEachSideOfBorder(permuted$peeps1, neighbours, neighbours.queen, permuted.sp)
  
  #Multiply those two
  mean(a * b$values)
  
})

# }

proc.time()-x

plot(density(null.moran.neigh.min))


a.edges <- selectiveNeighbourIndex(gridedges$peeps1, neighbours, useNeighMax = F)
b.edges <- localMoransI_onEachSideOfBorder(gridedges$peeps1, neighbours, neighbours.queen, gridedges.sp)

#Multiply those two
edges.moran.neigh.min <- mean(a.edges * b.edges$values)


a.smooth <- selectiveNeighbourIndex(gridsmooth$peeps1, neighbours, useNeighMax = F)
b.smooth <- localMoransI_onEachSideOfBorder(gridsmooth$peeps1, neighbours, neighbours.queen, gridsmooth.sp)

#Multiply those two
smooth.moran.neigh.min <- mean(a.smooth * b.smooth$values)


ggplot(data.frame(null = null.moran.neigh.min), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edges.moran.neigh.min, colour='red') +
  geom_vline(xintercept = smooth.moran.neigh.min, colour='green')




#~~~~~~~~~~~~~~~~~~~~~~
#AV IN VICINITY OF CELL * AACD----
#~~~~~~~~~~~~~~~~~~~~~~

#i.e. av welsh people in this cell and surrounding eight.
#i.e. for each border, it's gonna be i not j. So values will repeat for those borders.
#How to get those repeated for each border pair easily?

#Currently with standard 7*8 grids
#So this one is just plain neighbour av, which is:
#neighbours.queen <- poly2nb(gridsmooth.sp,queen = T)
#neighbours <- poly2nb(gridsmooth.sp,queen = F)
m <-  nb2mat(neighbours.queen)
#Actually, that's just repeating values a certain number of times. Must be easier way?
#This works I think...
#https://stackoverflow.com/questions/43186235/replicate-certain-values-in-vector-determined-by-other-vector

#i is repeated for each border. Sum should be 97...ah, no, it's 194 of course!
# repnumber <- sapply(neighbours,length)
# repz <- rep()

#Get null
x <- proc.time()
null.vicinity.aacd <- sapply(1:300, function(x) {
  
  cat(x,'\n')

  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  
  a <- avInVicinityForEachBorderPair(m,permuted$peeps1)
  b <- aacd(data.sf = permuted, var = 'peeps1')
  
  #Multiply those two
  mean(a * b$values)
  
})
proc.time()-x

plot(density(null.vicinity.aacd))


#SMOOTH
lag4calc <- avInVicinityForEachBorderPair(m,gridsmooth$peeps1)
#Get plain ol' ACD
acd <- aacd(data.sf = gridsmooth, var = 'peeps1')

#And final is just two multiplied / divided by total cell no
smooth.vicinity.aacd <- mean(acd$values * lag4calc)

#GRID
lag4calc <- avInVicinityForEachBorderPair(m,gridedges$peeps1)
#Get plain ol' ACD
acd <- aacd(data.sf = gridedges, var = 'peeps1')

#And final is just two multiplied / divided by total cell no
edge.vicinity.aacd <- mean(acd$values * lag4calc)


ggplot(data.frame(null = null.vicinity.aacd), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edge.vicinity.aacd, colour='red') +
  geom_vline(xintercept = smooth.vicinity.aacd, colour='green')


data.frame(acd = acd$values, lag = lag4calc) %>% 
  gather(key = thing, value = value) %>% 
  ggplot(aes(x = thing, y = value)) +
  geom_boxplot()


#Look at both grid types one one plot
# gridboth <- gridedges
# gridboth <- gridboth %>% rename(peeps1edges = peeps1) %>% dplyr::select(-peeps2)
# gridboth$peeps1smooth <- gridsmooth$peeps1
# plot(gridboth)



#~~~~~~~~~~~~~~~~~~~~~~
#AV IN VICINITY OF CELL * SELECTIVE NEIGHBOUR (MAX)----
#~~~~~~~~~~~~~~~~~~~~~~

#i.e. av welsh people in this cell and surrounding eight.
#i.e. for each border, it's gonna be i not j. So values will repeat for those borders.
#How to get those repeated for each border pair easily?

#Currently with standard 7*8 grids
#So this one is just plain neighbour av, which is:
#neighbours.queen <- poly2nb(gridsmooth.sp,queen = T)
#neighbours <- poly2nb(gridsmooth.sp,queen = F)
m <-  nb2mat(neighbours.queen)
#Actually, that's just repeating values a certain number of times. Must be easier way?
#This works I think...
#https://stackoverflow.com/questions/43186235/replicate-certain-values-in-vector-determined-by-other-vector

#i is repeated for each border. Sum should be 97...ah, no, it's 194 of course!
# repnumber <- sapply(neighbours,length)
# repz <- rep()

#Get null
x <- proc.time()
null.vicinity.neighmax <- sapply(1:1000, function(x) {
  
  cat(x,'\n')
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  
  a <- avInVicinityForEachBorderPair(m,permuted$peeps1)
  b <- selectiveNeighbourIndex(permuted$peeps1,neighbours, useNeighMax = T)
  
  #Multiply those two
  mean(a * b)
  
})
proc.time()-x

plot(density(null.vicinity.aacd))


#SMOOTH
a <- avInVicinityForEachBorderPair(m,gridsmooth$peeps1)
b <- selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, useNeighMax = T)

#And final is just two multiplied / divided by total cell no
smooth.vicinity.neighmax <- mean(a*b)

#GRID
a <- avInVicinityForEachBorderPair(m,gridedges$peeps1)
b <- selectiveNeighbourIndex(gridedges$peeps1,neighbours, useNeighMax = T)

#And final is just two multiplied / divided by total cell no
edge.vicinity.neighmax <- mean(a*b)


ggplot(data.frame(null = null.vicinity.neighmax), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edge.vicinity.neighmax, colour='red') +
  geom_vline(xintercept = smooth.vicinity.neighmax, colour='green')







#~~~~~~~~~~~~~~~~~~~~~~
#AV IN VICINITY OF CELL * RCI * AACD----
#~~~~~~~~~~~~~~~~~~~~~~

m <-  nb2mat(neighbours.queen)

#Get null
x <- proc.time()
null.vicinity.rci.aacd <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  
  a <- avInVicinityForEachBorderPair(m,permuted$peeps1)
  b <- aacd(data.sf = permuted, var = 'peeps1')
  c <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "peeps1")
  
  #Multiply those two
  mean(a * b$values * c$values)
  
})
proc.time()-x

plot(density(null.vicinity.rci.aacd))


#SMOOTH
a <- avInVicinityForEachBorderPair(m,gridsmooth$peeps1)
b <- aacd(data.sf = gridsmooth, var = 'peeps1')
c <- centralisationIndex_onEachSideOfBorder(data.sf = gridsmooth, var = "peeps1")

#And final is just two multiplied / divided by total cell no
smooth.vicinity.rci.aacd <- mean(a * b$values * c$values)

#GRID
a <- avInVicinityForEachBorderPair(m,gridedges$peeps1)
b <- aacd(data.sf = gridedges, var = 'peeps1')
c <- centralisationIndex_onEachSideOfBorder(data.sf = gridedges, var = "peeps1")

#And final is just two multiplied / divided by total cell no
edge.vicinity.rci.aacd <- mean(a * b$values * c$values)



ggplot(data.frame(null = null.vicinity.rci.aacd), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edge.vicinity.rci.aacd, colour='red') +
  geom_vline(xintercept = smooth.vicinity.rci.aacd, colour='green')





#~~~~~~~~~~~~~~~~~~~~~~
#AV IN VICINITY OF CELL * RCI * SELECTIVE NEIGHBOUR (MAX)----
#~~~~~~~~~~~~~~~~~~~~~~

m <-  nb2mat(neighbours.queen)

#Get null
x <- proc.time()
null.vicinity.rci.neighmax <- sapply(1:100, function(x) {
  
  cat(x,'\n')
  
  permuted <- gridedges
  permuted$peeps1 <- sample(permuted$peeps1,nrow(permuted),replace = F)
  
  a <- avInVicinityForEachBorderPair(m,permuted$peeps1)
  b <- selectiveNeighbourIndex(permuted$peeps1,neighbours, useNeighMax = T)
  c <- centralisationIndex_onEachSideOfBorder(data.sf = permuted, var = "peeps1")
  
  #Multiply those two
  mean(a * b * c$values)
  
})
proc.time()-x

plot(density(null.vicinity.rci.neighmax))


#SMOOTH
a <- avInVicinityForEachBorderPair(m,gridsmooth$peeps1)
b <- selectiveNeighbourIndex(gridsmooth$peeps1,neighbours, useNeighMax = T)
c <- centralisationIndex_onEachSideOfBorder(data.sf = gridsmooth, var = "peeps1")

#And final is just two multiplied / divided by total cell no
smooth.vicinity.rci.neighmax <- mean(a * b * c$values)

#GRID
a <- avInVicinityForEachBorderPair(m,gridedges$peeps1)
b <- selectiveNeighbourIndex(gridedges$peeps1,neighbours, useNeighMax = T)
c <- centralisationIndex_onEachSideOfBorder(data.sf = gridedges, var = "peeps1")

#And final is just two multiplied / divided by total cell no
edge.vicinity.rci.neighmax <- mean(a * b * c$values)



ggplot(data.frame(null = null.vicinity.rci.neighmax), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = edge.vicinity.rci.neighmax, colour='red') +
  geom_vline(xintercept = smooth.vicinity.rci.neighmax, colour='green')











#~~~~~~~~~~~~~~~~~~~~~~
#CREATE LARGER GRID----
#~~~~~~~~~~~~~~~~~~~~~~

#But with similar patterning. Previous is 7*8
#Have to make in QGIS to keep sizes square. Double the size? 14*16

gridsmooth <- st_read('saves/shapefiles/7by8gridSquare.shp')
gridedges <- st_read('saves/shapefiles/7by8gridSquare.shp')

#Of course, it would be oriented completely differently. Can I use ID to fix?
plot(gridsmooth)

#Reorient via ID so order is top-left across
neworder <- list()
for(i in 1:7) neworder[[length(neworder)+1]] <- seq(from=i, to=56, by=7)

gridsmooth$id <- unlist(neworder)
gridsmooth <- gridsmooth %>% arrange(id)

gridedges$id <- unlist(neworder)
gridedges <- gridedges %>% arrange(id)



#Coding values directly from G (left to right)
#G uses single value - for DI for two pops, the other peeps are 1-x
gridsmooth$peeps1 <- c(0.05,0.2,0.4,0.6,0.8,0.8,0.8,
                       0.05,0.2,0.4,0.6,0.8,0.8,0.8,
                       0.05,0.2,0.4,0.6,0.6,0.8,0.8,
                       0.05,0.2,0.4,0.6,0.6,0.6,0.6,
                       0.05,0.2,0.4,0.4,0.4,0.4,0.4,
                       0.05,0.2,0.2,0.2,0.2,0.2,0.2,
                       0.05,0.2,0.2,0.2,0.2,0.2,0.2,
                       0.05,0.05,0.05,0.05,0.05,0.05,0.05)
gridsmooth$peeps2 <- 1 - gridsmooth$peeps1

gridedges$peeps1 <- c(0.2,0.4,0.4,0.4,0.4,0.4,0.2,
                      0.2,0.2,0.2,0.05,0.2,0.2,0.4,
                      0.2,0.2,0.05,0.8,0.05,0.6,0.4,
                      0.2,0.05,0.8,0.8,0.8,0.05,0.4,
                      0.2,0.05,0.8,0.8,0.8,0.05,0.4,
                      0.2,0.05,0.05,0.8,0.05,0.05,0.6,
                      0.2,0.6,0.05,0.05,0.05,0.6,0.2,
                      0.2,0.2,0.6,0.6,0.6,0.6,0.2)
gridedges$peeps2 <- 1 - gridedges$peeps1

plot(gridsmooth)
plot(gridedges)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PREPPING 16X16 SIM CITIES----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#First prepped in excel
#CSV contains three 16x16 cities stacked on top of one another:
#1. Smooth, from top down. 
#2 Edge with centre
#3. Edge, smooth cut in half and shifted 4 cells so 0 is next to 1
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


ncol=16
nrow=16
cutoff=0

#Check on optimised version for basic AACD.
#Orig version should still work for optimisation even if actual val is wrong (is double-counting)
smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
edgecentreAACD <- getAverageAbsoluteContiguousDifference(gridedge_centre$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
edgesplitAACD <- getAverageAbsoluteContiguousDifference(gridedge_split$peeps1, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)

actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD), cutoff = cutoff)

#will be same for any of the 3, same values...
permute16 = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff, torus=F)

#Err. I wanted to optimise didn't I?
ggplot(data.frame(null = permute16), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = smoothAACD, colour='red') +
  geom_vline(xintercept = edgecentreAACD, colour='green') +
  geom_vline(xintercept = edgesplitAACD, colour='blue')


#ALSO TRYING WITH SQUARED / OTHER EXPONENT VALUES TO SEE IF MAKES DIFF...newp! Well some, but...
expz=9

smoothAACD <- getAverageAbsoluteContiguousDifference(gridsmooth$peeps1^expz, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
edgecentreAACD <- getAverageAbsoluteContiguousDifference(gridedge_centre$peeps1^expz, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)
edgesplitAACD <- getAverageAbsoluteContiguousDifference(gridedge_split$peeps1^expz, ncol = ncol, nrow = nrow, cutoff = cutoff, torus = F)

# actualAACDs <- data.frame(type = c('smooth','edge'), AACD = c(smoothAACD,edgeAACD), cutoff = cutoff)

#will be same for any of the 3, same values...
permute16 = getRepeatedAACDfromPermutedCells(attribute = gridsmooth$peeps1^expz, ncol = ncol, nrow = nrow, numreps = 200000, cutoff = cutoff, torus=F)

#Err. I wanted to optimise didn't I?
ggplot(data.frame(null = permute16), aes(x=null)) +
  geom_density() +
  geom_vline(xintercept = smoothAACD, colour='red') +
  geom_vline(xintercept = edgecentreAACD, colour='green') +
  geom_vline(xintercept = edgesplitAACD, colour='blue')





#Min and max...
min <- optimiseAverageAbsoluteContiguousDifference(attribute = gridsmooth$peeps1, 
                                                   secondpop = 1-gridsmooth$peeps1,
                                                   ncol = ncol, nrow = nrow, cutoff = 0, maximise = F, breakval = 1000000)


#Min from random - will be different optimised pattern. Smooth is already min, nowhere to go...
gridsmooth$randomised <- gridsmooth$peeps1[sample(1:nrow(gridsmooth),replace=F)]
minfromrandom <- optimiseAverageAbsoluteContiguousDifference(attribute = gridsmooth$randomised, 
                                                   secondpop = 1-gridsmooth$randomised,
                                                   ncol = ncol, nrow = nrow, cutoff = 0, maximise = F, breakval = 3000000)

max <- optimiseAverageAbsoluteContiguousDifference(attribute = gridsmooth$peeps1, 
                                                   secondpop = 1-gridsmooth$peeps1,
                                                   ncol = ncol, nrow = nrow, cutoff = 0, maximise = T, breakval = 1000000)

minaacd = getAverageAbsoluteContiguousDifference(min$attribute, ncol = ncol, nrow = nrow, cutoff=0)
maxaacd = getAverageAbsoluteContiguousDifference(max$attribute, ncol = ncol, nrow = nrow, cutoff=0)

gridopt <- gridsmooth
gridopt$min <- min$attribute
gridopt$minfromrandom <- minfromrandom$attribute
gridopt$max <- max$attribute

plot(gridopt[,'min'])#Makes sense: this is likely the absolute minimum.
plot(gridopt[,'minfromrandom'])
plot(gridopt[,'max'])#Yup!


#~~~~~~~~~~~~~~
#VARIOGRAM?----
#~~~~~~~~~~~~~~

#https://stats.idre.ucla.edu/r/faq/how-do-i-generate-a-variogram-for-spatial-data-in-r/
library(geoR)

#http://www.leg.ufpr.br/geoR/geoRdoc/geoRintro.html
data(s100)

#OK, so how to make our smooth data into something this can read?
smooth.geor <- gridsmooth %>% st_set_geometry(NULL)
st_coordinates(gridsmooth)[,c(1,2)]

#Add coordinates (should be centroids, I think...)
smooth.geor <- cbind(gridsmooth %>% st_set_geometry(NULL),st_coordinates(st_centroid(gridsmooth))[,c(1,2)])

edge1.geor <- cbind(gridedge_centre %>% st_set_geometry(NULL),st_coordinates(st_centroid(gridedge_centre))[,c(1,2)]) 

edge2.geor <- cbind(gridedge_split %>% st_set_geometry(NULL),st_coordinates(st_centroid(gridedge_split))[,c(1,2)]) 

#chessboard...
chessboard.geor <- cbind(gridopt %>% st_set_geometry(NULL) %>% dplyr::select(max),
                         st_coordinates(st_centroid(gridopt))[,c(1,2)]) 


#max distance for variogram is...
st_bbox(gridsmooth)[3]-st_bbox(gridsmooth)[1]

#Now, um...
x <- variog(coords = smooth.geor[,c('X','Y')], data = smooth.geor$peeps1, uvec=seq(0,8000,l=50))
plot(x)

y <- variog(coords = edge1.geor[,c('X','Y')], data = edge1.geor$peeps1, uvec=seq(0,8000,l=50))
plot(y)

a <- variog(coords = edge2.geor[,c('X','Y')], data = edge2.geor$peeps1, uvec=seq(0,8000,l=50))
plot(a)

z <- variog(coords = chessboard.geor[,c('X','Y')], data = chessboard.geor$max, uvec=seq(0,8000,l=50))
plot(z)

#Oh, needs to be ggplot. Harrumph.
#cowplot::plot_grid(plotlist = list(x,y,a,z), labels = c('smooth','edge centre','edge split','chessboard'))






#Looking at the full data... Well, not much use
x <- variog(coords = smooth.geor[,c('X','Y')], option='cloud',data = smooth.geor$peeps1, uvec=seq(0,8000,l=50))
plot(x)

y <- variog(coords = edge1.geor[,c('X','Y')], option='cloud', data = edge1.geor$peeps1, uvec=seq(0,8000,l=50))
plot(y)

a <- variog(coords = edge2.geor[,c('X','Y')], option='cloud', data = edge2.geor$peeps1, uvec=seq(0,8000,l=50))
plot(a)

z <- variog(coords = chessboard.geor[,c('X','Y')], option='cloud', data = chessboard.geor$max, uvec=seq(0,8000,l=50))
plot(z)


#Boxplots might be better... hmm, nah?
x <- variog(coords = smooth.geor[,c('X','Y')], data = smooth.geor$peeps1, uvec=seq(0,8000,l=50), bin.cloud=T)
plot(x)
plot(x, bin.cloud=T)

z <- variog(coords = chessboard.geor[,c('X','Y')], data = chessboard.geor$max, uvec=seq(0,8000,l=50), bin.cloud=T)
plot(z, bin.cloud=T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#MAKE OWN VARIOGRAM WITH ABSOLUTE DIFFS----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Using the same cities as before. For this, we're gonna need a full distance matrix
#And a way to convert to pairs to then bin distances to find absolute differences
#Something I did somewhere a while back, but let's start from scratch...

#Shifted all to function
x=acd_variogram(gridsmooth,'peeps1')
y=acd_variogram(gridedge_centre,'peeps1')
z=acd_variogram(gridedge_split,'peeps1')
#Chessboard pattern
a=acd_variogram(gridopt,'max')
#Minimised from random (so new smooth pattern)
b=acd_variogram(gridopt,'minfromrandom')

plot_grid(plotlist = list(x,y,z,a,b),labels=c('smooth','edge centre','edge split','chessboard','min from random'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ADD EXTRA LAG TO CENTRALISATION INDEX----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#testing how I think that works
gridsmooth.sp <- as_Spatial(gridsmooth)
gridsmooth.neighbours <- poly2nb(gridsmooth.sp, queen=T)

#https://rdrr.io/rforge/spdep/man/nblag.html
gridsmooth.neighbours.2 <- nblag(gridsmooth.neighbours,maxlag = 3)
gridsmooth.neighbours.2 <- nblag_cumul(gridsmooth.neighbours.2)#Combine into one list

#Those two are same
gridsmooth.neighbours.2[[25]]
gridsmooth.neighbours[[25]]

#Second index contains 2nd order lag
gridsmooth.neighbours.2[[2]]

#So in theory, to get specific set of neighbours and lag, collapse both?
c(gridsmooth.neighbours.2[[1]][[25]],gridsmooth.neighbours.2[[2]][[25]])

#Is that right? Yes, although... how is that 25?? Is it upside down??
gridsmooth$check <- 0
gridsmooth$check[gridsmooth.neighbours.2[[150]]] <- 1

#Is that right? Yes, although... how is that 25?? Is it upside down?? Would appear so. Huh.
gridsmooth$check <- 0
gridsmooth$check[25] <- 1
gridsmooth$check[gridsmooth$id==25] <- 1

plot(gridsmooth)

#OK, that's the code anyway. Now to add to centralisation thingyo.









