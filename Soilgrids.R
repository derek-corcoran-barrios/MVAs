library(rgee)
ee_Initialize()
db <- "projects/soilgrids-isric/ocs_mean"
image <- ee$Image(db)$select("ocs_0-30cm_mean")

geometry <- ee$Geometry$Rectangle(
  coords = c(-74.84180, -44.08860, -71.58262, -39.35440),
  proj = "EPSG:4326",
  geodesic = FALSE
)

a <- ee_as_raster(image = image,
                  region = geometry,
                  via = "drive")

a <- projectRaster(a, crs = "+proj=longlat +datum=WGS84 +no_defs")

saveRDS(a, "NewCarbonStock.rds")

a <- a %>% crop(Regiones)
a <- a %>% mask(rasterize(Regiones, a))

values(a) <- ifelse(values(a) == 0, NA, values(a))
