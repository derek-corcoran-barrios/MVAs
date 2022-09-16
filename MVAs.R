library(tidyverse)
library(raster)
library(readxl)
      

  
DF <- read_excel("PaleoNicheModeling/BodyMassSDM2020.xlsx", 
                   sheet = "Summary") %>% 
        rename(Mass = `Body mass (kg)`, Feeding = Diet) %>% 
        mutate(Minimum = 3876, Minimum_Lo = 5095) %>%
        # Herbivores = 2.71*(Mass^1.02)*0.01, 
        # Omnivores = 3.4*(Mass^0.92)*0.01,
        # Carniivores = 137*(Mass^1.37)*0.01,
   
        mutate(HomeRange = case_when(Feeding == "Herbivore" ~ 2.71*(Mass^1.02)*0.01,
                                     Feeding == "Omnivore" ~ 3.4*(Mass^0.92)*0.01,
                                     Feeding == "Carnivore" ~ 137*(Mass^1.37)*0.01),
               MinimumArea = Minimum*HomeRange, MinimumAreaLo = Minimum_Lo*HomeRange, Species = str_replace_all(Species, pattern = " ", replacement = "_")) %>% 
        dplyr::select(-Minimum, -Minimum_Lo, -Reference) %>% 
        arrange(Species) %>% 
    mutate(Species = str_replace_all(Species, "Odocoileus_virginanus", "Odocoileus_virginianus")) %>%
  mutate(Species = str_replace_all(Species,"Tapirus_pichanque", "Tapirus_pinchaque")) %>% 
  mutate(Species = str_replace_all(Species,"Tapirus_bairdii", "Tapirus_bairdi"))
  
      
### Area minima segun TSS
      
Files <- data.frame(File = list.files(path = "PaleoNicheModeling/PredictedTSS",full.names = T)) %>% 
  mutate(Species = str_remove_all(File, "PaleoNicheModeling/PredictedTSS/"), Species = str_remove_all(Species, "_TSS.rds"), Species = str_remove_all(Species, "\\("), Species = str_remove_all(Species, "\\)"), 
         Species = str_remove_all(Species, "_=Scelidodon_chiliense"),
         Species = str_replace_all(Species, "Holmesina_paulacoutoi", "Holmesina_paulacouti"),
         Species = str_replace_all(Species, "Pampatherium_humboldtii", "Pampatherium_humboldti"),
         File = as.character(File))

DF$Species[!(DF$Species %in% Files$Species)]

DF <- DF %>% left_join(Files) %>% dplyr::filter(!is.na(File))

Areas <- list()
      
AreaTemp <- area(read_rds(DF$File[1]))
      
  for(i in 1:nrow(DF)){
    Areas[[i]] <- data.frame(Species = DF$Species[i], year = seq(21000, 0, by = -1000), Distribution = NA)
        
    # area TSS
        
    RasterTemp <- read_rds(DF$File[i])
    Areas[[i]]$Distribution <- (AreaTemp*RasterTemp) %>% cellStats("sum")
        
    ## Area low Tri
        
    message(paste(i, "of", nrow(DF)))
    gc()
}
      
  Areas <- Areas %>% reduce(bind_rows) %>% full_join(DF) %>% mutate(NVA = Distribution/MinimumArea, NVA_Lo = Distribution/MinimumAreaLo) 
      
saveRDS(Areas, "Areas.rds")
write_csv(Areas, "Areas.csv")
### Clustering por años 

For_Clust <- Areas %>% 
  dplyr::select(Taxon_name, year, NVA) %>% 
  pivot_wider(names_from = year, values_from = NVA) 

Spp <- For_Clust$Taxon_name

For_Clust <- For_Clust %>% dplyr::select(-Taxon_name) %>% as.matrix()

rownames(For_Clust) <- Spp
scaled_data = as.matrix(scale(For_Clust))

cl <- kmeans(scaled_data, 4)

Grupos <-  data.frame(Group = cl$cluster, Taxon_name = Spp)

Areas_G <- full_join(Areas, Grupos)  %>% 
  group_split(Group) %>% 
  walk(~print(ggplot(.x, aes(x = year, y = NVA, group = Taxon_name)) + geom_ribbon(aes(fill = Taxon_name, ymax = NVA, ymin = NVA_Lo), alpha = 0.5) + geom_line(aes(color = Taxon_name)) + theme_bw() + theme(legend.position = "none") + geom_hline(yintercept = 1, lty = 2)+ facet_wrap(~Taxon_name)))


Extinct_Spp <- Areas %>% dplyr::filter(NVA < 1)
## Select number of groups

library(mclust)

  d_clust <- Mclust(as.matrix(scaled_data), G=1:15, 
                  modelNames = mclust.options("emModelNames"))
d_clust$BIC
plot(d_clust)

require("TSclust")
proxy::pr_DB$set_entry(FUN = diss.ACF, names = c("ACFD"),loop = TRUE, type = "metric", distance = TRUE,description = "Autocorrelation-based distance")



require("TSclust")
proxy::pr_DB$set_entry(FUN = diss.ACF, names = c("ACFD"),loop = TRUE, type = "metric", distance = TRUE,description = "Autocorrelation-based distance")

DTWndtw <- function(x, y, ...) {dtw(x, y, ...,step.pattern = asymmetric,distance.only = TRUE)$normalizedDistance}# Register the distance with proxy
proxy::pr_DB$set_entry(FUN = ndtw, names = c("nDTW"),loop = TRUE, type = "metric", distance = TRUE,description = "Normalized, asymmetric DTW")# Partitional clusteringtsclust(CharTraj[1L:10L], k = 2L,distance = "nDTW", seed = 838)