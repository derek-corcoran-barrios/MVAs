setwd("/home/derek/Documents/PaleoPlants")
rm(list=ls())

require("ncdf4")
require("lattice")
require("ggplot2")
library(raster)
library(tidyverse)
library(rworldxtra)
library(sf)
library(terra)
data("countriesHigh")

e <- readRDS("/home/derek/Documents/MVAs/PaleoNicheModeling/DiversityEstimates/Extinct only/HQ_Extinct.rds") %>% extent()

World <- countriesHigh %>% st_as_sf(countriesHigh) %>% st_make_valid() %>% st_crop(e)

file <- "LateQuaternary_Environment.nc";

env_nc      <- ncdf4::nc_open(file)
longitude   <- ncdf4::ncvar_get(env_nc, "longitude")
latitude    <- ncdf4::ncvar_get(env_nc, "latitude")
years       <- ncdf4::ncvar_get(env_nc, "time")
biome       <- ncdf4::ncvar_get(env_nc, "biome")

ncdf4::nc_close(env_nc)

#Table_Biomes <- tribble(
#                         ~Code,  ~Biome,  ~Color,
#                          0, "Water bodies", "#7fb3ff",
#                          1, "Tropical evergreen forest", "#64ff01"
#                          2,
#)

Table_Biomes <- data.frame(Code = c(0:28),
                           Biome = c("Water bodies",
                                     "Tropical evergreen forest",
                                     "Tropical semi-deciduous forest",
                                     "Tropical deciduous forest/woodland",
                                     "Temperate deciduous forest",
                                     "Temperate conifer forest",
                                     "Warm mixed forest",
                                     "Cool mixed forest",
                                     "Cool conifer forest",
                                     "Cold mixed forest",
                                     "Evegreen taiga/montane forest",
                                     "Deciduous taiga/montane forest",
                                     "Tropical savanna",
                                     "Tropical xerophytic shrubland",
                                     "Temperate xerophytic shrubland",
                                     "Temperate sclerophyll woodland",
                                     "Temperate broadleaved savanna",
                                     "Open conifer woodland",
                                     "Boreal parkland",
                                     "Tropical grassland",
                                     "Temperate grassland",
                                     "Desert",
                                     "Steppe tundra",
                                     "Shrub tundra",
                                     "Dwarf shrub tundra",
                                     "Prostrate shrub tundra",
                                     "Cushion forb lichen moss tundra",
                                     "Barren",
                                     "Land ice"))

setwd("/home/derek/Documents/MVAs/")

model.matrix()

library(gifski)
Tabla <- list() 

Tiempos <- seq(-21000, 0, by = 1000)

Pallete <- c("Barren" = "#b9b5a9", "Cool conifer forest" = "#00a20d", "Cool mixed forest" = "#caffc8", "Cushion forb lichen moss tundra" = "#9368eb", 
             "Deciduous taiga/montane forest" = "#67aaf8", "Desert" = "#f3fcc5", "Dwarf shrub tundra" = "#7a843f", 
             "Evegreen taiga/montane forest" = "#0523b7", "Land ice" = "#b4cfda", "Open conifer woodland" = "#ff97e9", 
             "Prostrate shrub tundra" = "#d3b0aa", "Shrub tundra" = "#5cff8c", "Steppe tundra" = "#f2e92a", "Temperate conifer forest" = "#19856e", 
             "Temperate grassland" = "#ffebab", "Temperate sclerophyll woodland" = "#8b982e", "Temperate xerophytic shrubland" = "#fbe0cb", 
             "Tropical deciduous forest/woodland" = "#ad802b", "Tropical evergreen forest" = "#184411", 
             "Tropical grassland" = "#f3b03d", "Tropical savanna" = "#b2ed45", "Tropical semi-deciduous forest" = "#629200", 
             "Tropical xerophytic shrubland" = "#fdb58f", "Warm mixed forest" = "#070761", "Water bodies" = "black")

save_gif(expr = for(i in 1:length(Tiempos)){
  Time <- biome[,,years == Tiempos[i]]
  Temp <- raster(t(Time)[300:1,], 
                 xmn = min(longitude), xmx = max(longitude),
                 ymn = min(latitude), ymx = max(latitude))
  
  Area <- Temp %>% 
    raster::area() %>% 
    crop(e) %>% 
    as("SpatialPixelsDataFrame") %>% 
    as.data.frame() %>% 
    rename(Area = layer)
  
  Temp <- Temp %>% 
    crop(e) %>% 
    as("SpatialPixelsDataFrame") %>% 
    as.data.frame() %>% rename(Code = layer) %>% 
    left_join(Table_Biomes) %>% 
    left_join(Area)
  
  Tabla[[i]] <- Temp %>% mutate(Year = -1*Tiempos[i])
  
  p <- ggplot() + geom_raster(data =Temp, aes(x = x, y = y, fill = Biome)) +
    geom_sf(data = World, alpha = 0) +
    theme_bw() + theme(legend.position = "bottom") +
    labs(x = NULL,
         y = NULL,
         title = paste(-1*Tiempos[i], "YBP")) + scale_fill_manual(name = "biome",values = Pallete)
  
  print(p)
}, "Biomes.gif")



Tabla <- Tabla %>% reduce(bind_rows)

Areas_Biomas <- Tabla %>% 
  group_by(Biome, Year) %>% 
  summarise(Area = sum(Area))

ggplot(data = Areas_Biomas, aes(x = Area, y = Year, fill = Biome)) + 
  geom_area() + 
  theme(legend.position = "none")

g <- ggplot() + geom_raster(data =Tabla, aes(x = x, y = y, fill = Biome)) +
  geom_sf(data = World, alpha = 0) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = NULL,
       y = NULL) + facet_wrap(~Year)


png(filename = "Biomas.png", width = 5000, height = 5000, res = 300)
print(g)
dev.off()


Tabla2 <- Tabla %>% dplyr::filter(Year %in% c(0,  2000, 4000,  6000, 8000,10000, 12000, 14000, 16000, 18000, 
                                              20000, 21000))



g <- ggplot() + geom_raster(data =Tabla2, aes(x = x, y = y, fill = Biome)) +
  geom_sf(data = World, alpha = 0) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = NULL,
       y = NULL) + facet_wrap(~Year) + scale_fill_manual(name = "biome",values = Pallete)


png(filename = "Biomas2.png", width = 3500, height = 3500, res = 300)
print(g)
dev.off()


Tabla2$Biome %>% unique() %>% as.character() %>% sort() %>% dput()

