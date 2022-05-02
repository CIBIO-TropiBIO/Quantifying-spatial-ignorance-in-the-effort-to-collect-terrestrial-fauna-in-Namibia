######################
# Mapping ignorance
# Author: Ricardo Correia
# Last updated: 15/07/2021
# Downloading records from GBIF for Angola and mapping ignorance
######################

# Load relevant libraries
library(rgbif)
library(rnaturalearth)
library(BIRDS) #Este paquete ya no está en CRAN, hay que descargarse la última version
##Hay que instalar las dependencias ‘dbscan’, ‘geosphere’, ‘mapedit’, ‘rgeos’, ‘shotGroups’, ‘taxize’
##install.packages(c('dbscan','geosphere','mapedit','rgeos','shotGroups','taxize'))
library(sf)

# Check how many georeferrenced bird records in Angola
#occ_count: es del paquete rgbif y te devuelve el numero de registros de ocurrencia
occ_count(taxonKey = 734,
          country = isocodes[grep("Angola", isocodes$name), "code"],
          georeferenced = T)

# Get records for birds of Angola
bats_angola <- occ_search(taxonKey = 734,
                          country = isocodes[grep("Angola", isocodes$name), "code"],
                          hasCoordinate=T,
                          fields = 'all',
                          limit = 100000) 

# Get country limits polygon
#Hay que instalar "rnaturalearthdata"
angola_poly <- ne_countries(scale = 'medium',
                            country = 'Angola',
                            returnclass = 'sf')
#angola_poly tiene que ser un sf
#ne_countries devuelve los poligonos de paises a una escala determinada
#Lo hacemos para Namibia
namibia_poly <- ne_countries(scale = 'medium',
                            country = 'Namibia',
                            returnclass = 'sf')

# Make grid from polygon
makeGrid3 <- function (poly, gridSize, hexGrid = TRUE, offset = NULL, buffer = FALSE, 
                       simplify = FALSE, crs.trans = "default", crs.return = "default", tol = 0.01) 
{
  gridSizeM <- gridSize * 1000
  if (!any(class(poly) %in% c("sfc", "sf", "SpatialPolygons", 
                              "SpatialPolygonsDataFrame"))) {
    stop("Entered polygon is not an sf, SpatialPolygon nor SpatialPolygonsDataFrame")
  }
  if (any(class(poly) %in% c("SpatialPolygons", "SpatialPolygonsDataFrame"))) {
    poly <- st_as_sf(poly)
  }
  if (is.na(st_crs(poly))) {
    stop("The polygon has no coordinate projection system (CRS) associated")
  }
  if(crs.trans == "default"){
    poly <- st_transform(poly, crs = st_crs(getUTMproj(poly)))  
  } else {
    poly <- st_transform(poly, crs = crs.trans)
  }
  dif <- as.numeric(abs(diff(matrix(st_bbox(poly), ncol = 2))))
  if (any(gridSizeM >= dif)) {
    stop("Grid cells must be smaller than the sampling area")
  }
  if (any(gridSizeM <= dif/500)) {
    message("Grid cells are too many (>=500), this may result in very long computation times")
  }
  if (simplify) {
    poly <- st_simplify(poly, dTolerance = tol)
  }
  if (buffer) {
    poly <- st_buffer(poly, dist = gridSizeM)
  }
  if (is.null(offset)) {
    offset <- st_bbox(poly)[c("xmin", "ymin")]
  }
  else {
    if (length(offset) != 2 || !all(is.integer(offset)) || 
        !is.numeric(offset)) 
      stop("Offset should be either NULL or numeric of length 2; lower left corner coordinates (x, y) of the grid")
  }
  grid <- st_make_grid(poly, cellsize = gridSizeM, square = !hexGrid, 
                       offset = offset, what = "polygons")
  if(crs.return == "default"){
    grid <- st_transform(grid, crs = st_crs(4326))
  } else {
    grid <- st_transform(grid, crs = crs.return)
  }
  
  return(grid)
}

#angola_grid es un sfc_POLYGON que es una colección de geometrias. Single_Feature_Collection
angola_grid <- makeGrid3(poly = angola_poly,
                        hexGrid = F,
                        crs.trans = "+proj=eck4",
                        crs.return = "+proj=eck4",
                        gridSize = 100)
str(angola_grid)

#para Namibia
namibia_grid <- makeGrid3(poly = namibia_poly,
                         hexGrid = F,
                         crs.trans = "+proj=eck4",
                         crs.return = "+proj=eck4",
                         gridSize = 100)
str(namibia_grid)



# Convert data to visit based format
angola_bat_obs <- organizeBirds(bats_angola$data, sppCol = "scientificName", simplifySppName = TRUE)

# Check visits lacking date information
angola_bat_obs$spdf <- angola_bat_obs$spdf[!apply(angola_bat_obs$spdf[,2:4], 1, function(x) {any(is.na(x))}),]

# Summarize visits by grid
angola_bat_summary <- summariseBirds(angola_bat_obs, grid = angola_grid)

# Calculate ignorance
angola_bat_ign <- exposeIgnorance(angola_bat_summary$spatial$nObs, h = 5)

# Plot ignorance
palBWR <- leaflet::colorNumeric(c("navyblue", "white","red"), c(0, 1), 
                                na.color = "transparent")
plot(angola_poly$geometry, col="grey90", border = "grey90", lwd=1)#Esto te dibuja Angola
plot(namibia_poly$geometry, col="grey90", border = "grey90", lwd=1)#Esto te dibuja Namibia

plot(angola_grid$geometry, col="grey90", border = "grey90", lwd=1)

plot(angola_bat_summary$spatial$geometry, col=palBWR(angola_bat_ign), border = 'grey', add=TRUE)#Esto te dibuja la grid
plot(angola_poly$geometry, col=NA, border = "black", lwd=1, add=TRUE)

legend("bottomleft", legend=c(seq(0, 1, length.out = 5), "NA"),
       col = c(palBWR(seq(0, 1, length.out = 5)), "grey90"),
       title = "Ignorance \nnObs, \nO0.5=5", pch = 15, bty="n")

