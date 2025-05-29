library(ggplot2)
library(sf)
library(leaflet)

apalachacola_sf <- st_read(data/S_USA.AdministrativeForest.shp)

st_layers(data/S_USA.AdministrativeForest.shp)
