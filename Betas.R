library(changeRangeR)
library(raster)
library(rasterVis)
library(rgdal)
library(Matrix.utils)
library(tidyverse)
if(Sys.info()["sysname"]== "Windows") library (parallelsugar)
mc.cores=4
### determine where you want outputs
summaryBaseDir= paste0(getwd(), "/Betas")
#summaryBaseDir='/Volumes/cm2/changeRangerDemos/trees190/selec'
if(!file.exists(summaryBaseDir)) dir.create(summaryBaseDir)



envGrid= readRDS("~/Documents/MVAs/PaleoNicheModeling/DiversityEstimates/Extinct only/HQ_Extinct.rds")[[1]]

mydir=paste0(getwd(),'/Test')

sumDirs=setupSummaryDirectories(summaryBaseDir, optionalSubDirs=c('funcDiv','betaDiv'))

saveRDS(sumDirs,file=paste0(sumDirs$sumBaseDir,'/dirList.rds'))

dir.create(path = paste0(getwd(), "/Temp"))
for (scn in allScen){
  # do this separately for each scenario, so that scenarios are the units for applying the workflow. have to link across scenarios by species names
  allSpeciesMaps=data.frame(rasterFiles=c(list.files(paste0(mydir),recursive=T, full.names=T)), stringsAsFactors=F)
 	allSpeciesMaps$sp.names=file_path_sans_ext(basename(allSpeciesMaps$rasterFiles))
 	allSpeciesMaps <- allSpeciesMaps %>% dplyr::filter(str_detect(rasterFiles,scn)) %>% mutate(sp.names = str_remove_all(sp.names, paste0("_TSS_", scn)))
 	sp.ind=speciesIndexTable(allSpeciesMaps,sumDirs)
 	cell.ind=cellIndexTable(envGrid,nCellChunks=5,sumDirs)
 	cellBySpeciesMatrices(sumDirs$cbsDir,
 	                      rasterFiles=allSpeciesMaps$rasterFiles,
 	                      spNames=allSpeciesMaps$spNames,
 	                      scenario=scn,
 	                      envGrid=envGrid,
 	                      sp.ind=sp.ind,
 	                      cell.ind=cell.ind,
 	                      nCellChunks=5, # number of chunks to split the envGrid into
 	                      mc.cores=mc.cores,
 	                      overwrite=T)
}
