library(changeRangeR)
library(raster)
library(rasterVis)
library(rgdal)
library(Matrix.utils)
library(tidyverse)
if(Sys.info()["sysname"]== "Windows") library (parallelsugar)
mc.cores=4

Spp <- list.files(path = "PaleoNicheModeling/PredictedTSS/", pattern = ".rds", full.names = F) %>% 
  str_remove_all(".rds")

Spp_Files <- list.files(path = "PaleoNicheModeling/PredictedTSS/", pattern = ".rds", full.names = T)

allScen= paste0(21:0, "kybp")

allScen <- ifelse(str_count(allScen) == 5, paste0(0, allScen), allScen)

m <- c(-Inf, 0.9, NA, 0.9, Inf, 1)
m <- matrix(m, ncol = 3, byrow = T)

for(i in 1:length(Spp_Files)){
  dir.create("Test")
  Temp <- read_rds(Spp_Files[i])
  Temp <- reclassify(Temp, m)
  for(j in 1:nlayers(Temp)){
    dir.create(paste0(getwd(),"/Test/", allScen[j]))
    writeRaster(Temp[[j]], filename = paste0(getwd(),"/Test/", allScen[j], "/", Spp[i], ".tif"))
  }
  message(paste(i, "of", length(Spp_Files)))
}


summaryBaseDir= paste0(getwd(), "/supertest")
if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)

#load environnent with the reprojected projection (typically the one you used for modeling). You only need the raster grid that , not the actual layer values
envGrid= readRDS("~/Documents/MVAs/PaleoNicheModeling/DiversityEstimates/Extinct only/HQ_Extinct.rds")[[1]]
# folder of binary range rasters
myDir=paste0(getwd(),'/Test')

# shapefiles for plotting. This one comes preinstallted
library(rworldxtra)
library(sf)

data("countriesHigh")

world.shp2= countriesHigh %>% st_as_sf() %>% st_make_valid()%>% st_crop(extent(envGrid)) %>% as_Spatial()

sumDirs=setupSummaryDirectories(summaryBaseDir, optionalSubDirs=c('funcDiv','phyloDiv'))
# its a good idea to save this in case you need to start an analysis in the middle
saveRDS(sumDirs,file=paste0(sumDirs$sumBaseDir,'/dirList.rds'))
str(sumDirs)

sumDirs <- readRDS("~/Documents/MVAs/supertest/dirList.rds")

allSpeciesMaps=tibble(rasterFiles=list.files(paste0(myDir,'/21kybp'),
                                             recursive=T, full.names=T)) %>% 
  mutate(sp.names= rasterFiles %>% basename %>% file_path_sans_ext) %>%
  separate(sp.names,into=c('g','sp','t','s'),sep='_') %>% select(-t,-s) %>%
  unite(sp.names,g,sp)
# species index table. columns: species name, integer index
sp.ind=speciesIndexTable(allSpeciesMaps,sumDirs)
# cell index Table. columns: long, lat, cellid

cellIndexTable2 <- function (env, nCellChunks, sumDirs, toInteger = T) 
{
  co = coordinates(env)
  keep = complete.cases(values(env))
  co = co[keep, ]
  if (toInteger) 
    co = apply(co, 2, as.integer)
  if (nCellChunks > 1) {
    chunks = cut(1:nrow(co), nCellChunks, labels = FALSE)
  }
  else {
    chunks = rep(1, nrow(co))
  }
  cell.ind = data.frame(co, cellID = as.integer(cellFromXY(env, 
                                                           co)), chunkID = as.integer(chunks))
  cell.ind =  cell.ind %>% data.frame
  saveRDS(cell.ind, file = paste0(sumDirs$sumBaseDir, "/cellIndexTable.rds"))
  cell.ind
}


cell.ind=cellIndexTable2(envGrid,nCellChunks=10,sumDirs, toInteger = F)

cell.ind=readRDS(paste0(sumDirs$myBaseDir,'/cellIndexTable.rds'))
head(cell.ind)



sp.ind=readRDS(paste0(sumDirs$myBaseDir,'/speciesIndexTable.rds'))
head(sp.ind)


cell.ind=readRDS(paste0(sumDirs$myBaseDir,'/cellIndexTable.rds'))
cell.ind = cell.ind[complete.cases(cell.ind),]

chunks.r=cellIDToRaster(cell.ind,envGrid,'chunkID')

fdMapPlot(chunks.r,shp=world.shp2,legend.args=list(text='chunk ID'))


for (scn in allScen){
  allSpeciesMaps=tibble(rasterFiles=list.files(paste0(myDir,'/',scn),
                                               recursive=T, full.names=T)) %>%
    mutate(spNames= rasterFiles %>% basename %>% file_path_sans_ext) %>%
      separate(spNames,into=c('g','sp','t','s'),sep='_') %>% select(-t,-s) %>%
    unite(spNames,g,sp)
  cellBySpeciesMatrices(sumDirs$cbsDir,
	                      rasterFiles=allSpeciesMaps$rasterFiles,
                        spNames=allSpeciesMaps$spNames,
	                      scenario=scn,
	                      envGrid=envGrid,
	                      sp.ind=sp.ind,
                        cell.ind=cell.ind,
	                      nCellChunks=10, # number of chunks to split the envGrid into
  	                      mc.cores=mc.cores,
	                      overwrite=T)
}
rich=lapply(allScen,function(scn){
  r=richnessFromCBS(cbsDir=sumDirs$cbsDir,
                    scenario=scn,env=envGrid,
                    mc.cores=mc.cores, outDir=sumDirs$richDir)
  # I like to store all the plots as I go
fdMapPlot(stack(r),paste0(sumDirs$figDir,'/Richness_',scn,'.pdf'),shp=world.shp2,
            legend.args=list(text='Species Richness'))
  r
}) %>% stack
fdMapPlot(rich,shp=world.shp2,legend.args=list(text='Species Richness'))


sp.ind=readRDS(paste0(sumDirs$sumBaseDir,'/speciesIndexTable.rds'))
ra=lapply(allScen,function(scn){
	    rangeArea(cbsDir=sumDirs$cbsDir,scenario=scn,
	              sp.ind=sp.ind,outDir=sumDirs$rangeSizeDir,mc.cores=mc.cores)
})
str(ra)
hist(ra[[1]]$rangeArea)


newSpInd= sp.ind %>% full_join(ra[[1]],by=c('species','index')) %>%
  rename(rangeAreaPresent=rangeArea) %>%
  full_join(ra[[2]],by=c('species','index')) %>%
  rename(rangeArea8580=rangeArea)
str(newSpInd)


### Rarity

#Here's an example using a 'species attribute', which is a common application. Here we use species attributes to refer to any property of an individual species. Common operations are to (1) summarize species attributes within a spatial unit, or (2) group species based on different values of an attribute and calculate a summary statistic. As an example of (1) we map species rarity (average value of 1/range area over species within a cell).
#Here's an example where we add a new place to store outputs, since rarity probably isn't common enough to include as a default.


sumDirs$rarityDir=file.path(sumDirs$rangeSizeDir,'Rarity')
if(!file.exists(sumDirs$rarityDir)) dir.create(sumDirs$rarityDir)
rar=lapply(allScen,function(scn){
  # generate the value of 1/range size for each species in an attribute table
  raritySpAttr=sumDirs$rangeSizeDir %>%
    list.files(full.names=T,pattern=scn) %>%
    readRDS %>%
    mutate(rarity=1/rangeArea) %>%
    select(-rangeArea)
  r=speciesAttributeByCell(cbsDir=sumDirs$cbsDir,scenario=scn,
                           attrTable=raritySpAttr, method='mean',
                           env=envGrid, outDir=sumDirs$rarityDir)
  fdMapPlot(log(r),plotFile=paste0(sumDirs$figDir,'/Rarity_',scn,'.pdf'),
            shp=world.shp2,legend.args=list(text='log(rarity)',line=2,side=4))
  r
}) %>% stack
fdMapPlot(log(rar),shp=world.shp2,legend.args=list(text='log(rarity)'))



compareRichness(sumDirs, allScen[1], allScen[2], plotFig = T) %>% plot()

for(i in 1:(length(allScen)-1)){
  Beta <- turnoverFromCBS(
    cbsDir = sumDirs$cbsDir,
    scn1 = allScen[i],
    scn2 = allScen[(i + 1)],
    env = envGrid,
    mc.cores = 2,
    betaDivChange = T,
    outputTable = F
  )
  saveRDS(Beta, paste0("Beta","_", allScen[i], "_", allScen[(i + 1)], ".rds"))
}

beepr::beep(8)

library(raster)
library(sf)
library(tidyverse)
library(tidymodels)


BetaNames <-  list.files(pattern = "Beta_") %>% str_remove_all("Beta_") %>% str_remove_all(".rds")

Betas <- list.files(pattern = "Beta_", full.names = T) %>% purrr::map(readRDS) %>% purrr::map(~.x[[7]]) %>% reduce(stack)


names(Betas) <- BetaNames

Betas_DF <- Betas %>% as("SpatialPixelsDataFrame") %>% as.data.frame()

colnames(Betas_DF)[1:21] <- BetaNames

Betas_DF <- Betas_DF %>% pivot_longer(contains("kybp"), names_to = "interval", values_to = "Beta_diversity") %>% dplyr::filter(!is.na(Beta_diversity))

Betas_DF <- Betas_DF %>% mutate(Interval_n = case_when(interval == "01kybp_00kybp" ~ 0.5,
                                                       interval == "02kybp_01kybp" ~ 1.5,
                                                       interval == "03kybp_02kybp" ~ 2.5,
                                                       interval == "04kybp_03kybp" ~ 3.5,
                                                       interval == "05kybp_04kybp" ~ 4.5,
                                                       interval == "06kybp_05kybp" ~ 5.5,
                                                       interval == "07kybp_06kybp" ~ 6.5,
                                                       interval == "08kybp_07kybp" ~ 7.5,
                                                       interval == "09kybp_08kybp" ~ 8.5,
                                                       interval == "10kybp_09kybp" ~ 9.5,
                                                       interval == "11kybp_10kybp" ~ 10.5,
                                                       interval == "12kybp_11kybp" ~ 11.5,
                                                       interval == "13kybp_12kybp" ~ 12.5,
                                                       interval == "14kybp_13kybp" ~ 13.5,
                                                       interval == "15kybp_14kybp" ~ 14.5,
                                                       interval == "16kybp_15kybp" ~ 15.5,
                                                       interval == "17kybp_16kybp" ~ 16.5,
                                                       interval == "18kybp_17kybp" ~ 17.5,
                                                       interval == "19kybp_18kybp" ~ 18.5,
                                                       interval == "20kybp_19kybp" ~ 19.5,
                                                       interval == "21kybp_20kybp" ~ 20.5))


library(mgcv)
Loess3d <- bam(Beta_diversity ~ s(Interval_n), data = Betas_DF)
# 21 -12- 8 - 0
# 21 -14- 7 - 0
# 21 -14- 7 - 0
# 21 - 16 - 10 -  5 -  0
# Conservador
# Relajado
# Relajado + Existentes
# Existentes

newDF <- expand.grid(Interval_n = seq(from = 21, to = 0, length.out = 200))


newDF$Pred <- predict(Loess3d, newDF)
newDF$SE <- predict(Loess3d, newDF, se.fit = T)$se.fit


ggplot(newDF, aes(x = Interval_n, y = Pred)) +
  geom_ribbon(aes(ymax = Pred + 3*SE, ymin = Pred - 3*SE), color = scales::muted("red"), alpha = 0.5) +
  geom_path(size = 0.5) +
  theme_bw() + labs(x ="kybp", y = "Beta diversity [Sørensen index]")

#### smoothing with tidymodels
library(tidymodels)

set.seed(2020)
data_split <- Betas_DF %>%
  initial_split(strata = Beta_diversity, prop = 0.5)
beta_train <- training(data_split)
beta_test  <- testing(data_split)

rm(data_split)
gc()

beta_rec <- 
  recipe(Beta_diversity ~ Interval_n, data = beta_train) %>% 
  step_ns(Interval_n, deg_free = tune("interval_df"))

beta_param <- 
  beta_rec %>% 
  parameters() %>% 
  update(
    "interval_df" = spline_degree()
  )

lm_mod <- linear_reg() %>% set_engine("lm")

set.seed(2020)
cv_splits <- vfold_cv(beta_train, v = 10, strata = Beta_diversity)

spline_grid <- grid_max_entropy(beta_param, size = 10)

beta_res <- tune_grid(lm_mod, beta_rec, resamples = cv_splits, grid = spline_grid)

saveRDS(beta_res, "betares.rds")
####

Betas_DF2 <- Betas_DF %>% dplyr::filter(interval %in% c("01kybp_00kybp",  "03kybp_02kybp",  
                                                        "05kybp_04kybp",  "07kybp_06kybp",  
                                                        "09kybp_08kybp",  "11kybp_10kybp",  
                                                        "13kybp_12kybp",  "15kybp_14kybp",  
                                                        "17kybp_16kybp",  "19kybp_18kybp", 
                                                        "21kybp_20kybp"))



g <- ggplot() + geom_raster(data =Betas_DF2, aes(x = x, y = y, fill = Beta_diversity)) +
  geom_sf(data = world.shp2, alpha = 0) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = NULL,
       y = NULL) + facet_wrap(~interval) + scale_fill_viridis_c(name = "Beta diversity")


png(filename = "betamap.png", width = 3000, height = 3000, res = 300)
print(g)
dev.off()

data("countriesHigh")

world.shp2= countriesHigh %>% st_as_sf() %>% st_make_valid()%>% st_crop(extent(envGrid))


ggplot(data = Betas_DF, aes(x = interval, y = Beta_diversity)) + 
  geom_boxplot() + 
  coord_flip() + 
  theme_bw()