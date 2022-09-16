library(raster)
library(sf)
library(tidyverse)
library(ggrepel)

## Extent de latitud y longitud

e <- new("Extent", xmin = -78.7573908010839, xmax = -57.2140532449429, 
         ymin = -58.5312247001755, ymax = -15.7199043140937)

Chile <- getData(name = "GADM", country = "CHL", level = 1) %>% 
  st_as_sf() %>% 
  st_crop(e)

Diversity <- read_rds("/home/derek/Documents/MVAs/PaleoNicheModeling/DiversityEstimates/Extinct only/HQ_Extinct.rds") %>%
  crop(Chile) 

Diversity2 <-Diversity %>% 
  mask(rasterize(Chile, Diversity)) %>% 
  as("SpatialPixelsDataFrame") %>%
  as.data.frame()

Diversity_DF <- Diversity2 %>% 
  pivot_longer(starts_with("Diversity"), names_to = "Year", values_to = "Richness") %>% 
  mutate(Year = as.numeric(str_remove_all(Year, "Diversity."))) %>% 
  dplyr::filter(Year %in% c(21000, 15000, 8000, 6000, 3000, 0))

ggplot() + 
  geom_raster(data = Diversity_DF, aes(x = x, y = y, fill = Richness)) + 
  geom_sf(data = Chile, alpha = 0, size = 0.2) + 
  theme_bw() +
  facet_wrap(~Year, nrow = 1) +
  scale_fill_viridis_c() +
  scale_x_continuous(breaks = c(-74, -70)) +
  labs(y = NULL,
       x = NULL)


### Supongamos que quieres hacer un zoom a alguna region:

## Elgies entre estas:

c("Aisén del General Carlos Ibáñez del Campo", "Antofagasta", 
  "Araucanía", "Arica y Parinacota", "Atacama", "Bío-Bío", "Coquimbo", 
  "Libertador General Bernardo O'Higgins", "Los Lagos", "Los Ríos", 
  "Magallanes y Antártica Chilena", "Maule", "Ñuble", "Región Metropolitana de Santiago", 
  "Tarapacá", "Valparaíso")

## Supongamos magallanes y aysen:

Patagonia <- Chile %>% dplyr::filter(NAME_1 %in% c("Magallanes y Antártica Chilena", "Aisén del General Carlos Ibáñez del Campo"))

Diversity <- read_rds("/home/derek/Documents/MVAs/PaleoNicheModeling/DiversityEstimates/Extinct only/HQ_Extinct.rds") %>%
  crop(Patagonia) 

Diversity2 <-Diversity %>% 
  mask(rasterize(Patagonia, Diversity)) %>% 
  as("SpatialPixelsDataFrame") %>%
  as.data.frame()

Diversity_DF <- Diversity2 %>% 
  pivot_longer(starts_with("Diversity"), names_to = "Year", values_to = "Richness") %>% 
  mutate(Year = as.numeric(str_remove_all(Year, "Diversity."))) %>% 
  dplyr::filter(Year %in% c(21000, 15000, 8000, 6000, 3000, 0))

ggplot() + 
  geom_raster(data = Diversity_DF, aes(x = x, y = y, fill = Richness)) + 
  geom_sf(data = Patagonia, alpha = 0, size = 0.2) + 
  theme_bw() +
  facet_wrap(~Year, nrow = 1) +
  scale_fill_viridis_c() +
  scale_x_continuous(breaks = c(-74, -70)) +
  labs(y = NULL,
       x = NULL)
    