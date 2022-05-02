library(sf)


# Get country limits polygon
# Install "rnaturalearthdata"
angola_poly <- ne_countries(scale = 'medium',
                            country = 'Angola',
                            returnclass = 'sf')



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
                         gridSize = 10)
str(angola_grid)

#para Namibia
namibia_grid <- makeGrid3(poly = namibia_poly,
                          hexGrid = F,
                          crs.trans = "+proj=eck4",
                          crs.return = "+proj=eck4",
                          gridSize = 10)
str(namibia_grid)

## Shapefiles

## Angola
## Directory where you want to save the shapefile
path = '/Users/javiermartinez/Documents/CIBIO/Proyectos/Eckert/Angola_shp/'
st_write(angola_grid, paste(path,'angola_grid_10km.shp'))


## Namibia
## Directory where you want to save the shapefile
path = '/Users/javiermartinez/Documents/CIBIO/Proyectos/Eckert/Namibia_shp/'
st_write(namibia_grid, paste(path,'namibia_grid_10km.shp'))
