library(dplyr)
library(tidyr)
library(sf)
library(rgbif)
library(viridis)
library(terra)
library(raster)
set.seed(1)

# Read in Shapefile for Study Area 
# Downloadable at BC maps at: https://bcgov.github.io/bcmaps/reference/bec.html
SA  <- sf::read_sf('F:/NewBap/SBS_SPS.shp')

### make a hexagon 
## Based on here: https://strimas.com/post/hexagonal-grids/
make_grid <- function(x, type, cell_width, cell_area, clip = FALSE) {
  if (!type %in% c("square", "hexagonal")) {
    stop("Type must be either 'square' or 'hexagonal'")
  }
  
  if (missing(cell_width)) {
    if (missing(cell_area)) {
      stop("Must provide cell_width or cell_area")
    } else {
      if (type == "square") {
        cell_width <- sqrt(cell_area)
      } else if (type == "hexagonal") {
        cell_width <- sqrt(2 * cell_area / sqrt(3))
      }
    }
  }
  # buffered extent of study area to define cells over
  ext1 <- vect(ext(x)+ cell_width) %>% as(., 'Spatial')
  crs(ext1) <- crs(x)
  # generate grid
  if (type == "square") {
    g <- rast(ext1, resolution = cell_width)
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    # generate array of hexagon centers
    g <- spsample(ext1, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    # convert center points to hexagons
    g_poly <- HexPoints2SpatialPolygons(g, dx = cell_width)
  }
  g_vect <-  vect(g_poly)
}

SA_hex <- make_grid(SA, type = 'hexagonal', cell_area = 100000000) # cell area is equal to 10,000 


## NEXT STEPS HERE
## 1. raster exctract most for the most column for each pixel 

grps <- rast('D:/Paper2_Clean/Spectral_Clusters/Structural_groups.tif')
library("sf")
library("exactextractr")

SA_hex_sf <- st_as_sf(SA_hex)  %>% st_set_crs(crs(SA))
###write a function to "hexify" data the spectral clusters  
hexify <- function(r, s){
  #step 1 calculcate the mode
  most_in_r <- exact_extract(r, s,  fun = 'mode')
  #s$mode <- most_in_r$mode
  mode_vect <- vect(s)
  s$ID = seq(1, nrow(s), by = 1)
  #step 2 calculate the coverages 
  r_NA <- r
  r_NA[!is.na(r)] <- 1
  #total coverage based on proportion that is not NA
  r_cov <- exact_extract(r_NA, s, fun = 'frac')
  df_r  <- data.frame(ID = seq(1, nrow (s)), 
                      mode = most_in_r, 
                      prop = r_cov)
  s_mod <- s %>% left_join(df_r)
  return(s_mod)
}
clust_hex  <- hexify(grps, SA_hex_sf)

## get the most common for each structural group 
most_in_A_alt <- exact_extract(grps, SA_hex_sf, fun = 'mode')
most_in_A <- terra::extract(grps, SA_hex_sf, fun ='modal', na.rm = T)
SA_hex_sf$mode <- most_in_A_alt
mode_vect <- vect(SA_hex_sf)

SA_hex_sf$ID <- seq(1, nrow(SA_hex_sf), by = 1)

# calculate the coverage for each structural group 
grps_NA <- grps
grps_NA[!is.na(grps_NA)] <- 1
grps_NA[is.na(grps_NA)] <-0

# Define the values to include in the vector
possible_values <- c("A", "B", "C")

# Generate a random vector of length 100 with values from possible_values
x <- sample(possible_values, 100, replace = TRUE)

## get frequency
get_fract <- function(x){
  counts <- table(x)
  sum_ <- sum(counts)
  most_common_prop <- (counts[which.max(counts)])/sum_
  return(most_common_prop)
  }


most_in_A_fract <- extract(grps_NA, SA_hex_sf, fun = get_fract) ### total coverved by a class
in_r <- rast('F:/Quesnel_RPA_Lidar2022/Clusters/SBS_clusters/Septemb_Clustsers/test_folder/structural_groups_clip.tif')
cropped <- crop(grps, SA_hex_sf) 
most_in_A_table <-exact_extract(grps, SA_hex_sf, c('variety', 'majority', 'frac'))
total_cover <- exact_extract(grps_NA, SA_hex_sf, c('frac'))
total_cover <- total_cover%>% 
most_in_A_table$ID = seq(1, nrow(most_in_A_table), by =1)
Extract_values <- most_in_A_table %>% filter(!is.na(majority)) %>% 
  pivot_longer(cols = starts_with('frac'),
               names_to = 'group', values_to = 'proportion') %>% 
  mutate(group = as.numeric(str_remove(group, pattern = 'frac_')),
         prop_max = ifelse(group == majority, proportion, NA)) %>% 
  group_by(ID) %>%
  dplyr::select(-c('proportion', 'group')) %>% distinct() %>% 
  group_by(ID) %>% 
  mutate(n_r = n()) %>% filter(!is.na(prop_max))

most_in_A_df <- data.frame(ID = most_in_A_table$ID) %>% 
  full_join(Extract_values, by = 'ID')
most_in_A_df$total_cover <- total_cover$frac_1
most_in_A_df$total_max <- most_in_A_df$total_cover * most_in_A_df$prop_max
SA_mode <- SA_hex_sf %>% cbind(most_in_A_df)
# SA_mode <- read_sf('F:/Quesnel_RPA_Lidar2022/Clusters/SBS_clusters/Septemb_Clustsers/Most_common.shp') %>% 
#   mutate(cov = as.numeric(cov)) %>% vect()
SA_mask <- vect(SA_mode) %>% mask(vect(SA)) 
clust_mask <- vect(clust_hex) %>% mask(vect(SA)) %>% st_as_sf()
SA_mode <- st_as_sf(SA_mask)

library(tmap)
tm_shape(SA_mode) + 
  tm_fill(col = 'mode', style = 'cat', showNA= F)

new_cols <- c("#a24936","#ffcf00","#777949","#006ba6","#E57128")

class(SA_mode)
write_sf(SA_mode, "D:/Paper2_Clean/Spectral_Clusters/Structural_Data/Most_common_nprop.shp")
library(maptiles)
SA_reprok <- project(SA_hex,"epsg:4326")
bg <- get_tiles(ext(SA_reprok), provider = 'Esri.WorldImagery', zoom = 7)
bg_reproj <- project(bg, crs(SA_hex)) %>% crop(SA_hex)
struct <- tm_shape(bg_reproj) + 
  tm_rgb() + 
  tm_shape(SA_mode) + 
  tm_fill(col = 'mode', palette = new_cols, showNA = F, 
          style = 'cat', colorNA = NULL,
          title = "Structural \nGroup")  + 
  tm_layout(legend.outside = F, 
            legend.outside.position = "left",
            legend.title.color = "white",
            legend.text.color = "white",
            legend.title.fontface = "bold",
            legend.position = c(0.8, 0)
  ) + 
  tm_scale_bar(position = c("LEFT", "TOP"))
## plot with DEM 
elev <- rast('F:/Sync/PhD_Writing/Paper2/QGIS_map/Spatial_AutoCorrMap/QGISLayers/elevation.tif')# %>% project(vect(SA_mode))

dem_map <- tm_shape(elev) + 
  tm_raster(palette = '-Greys', style = 'cont', 
            interpolate = T, midpoint = 900) +
  tm_shape(SA_mode) + 
  tm_fill(col = 'mode', palette = new_cols, showNA = F, 
          style = 'cat', colorNA = NULL,
          title = "Structural \nGroup")  + 
  tm_layout(legend.outside = F, 
            legend.outside.position = "left",
            legend.title.color = "white",
            legend.text.color = "white",
            legend.title.fontface = "bold",
            legend.position = c(0.8, 0)
  ) + 
  tm_scale_bar(position = c("LEFT", "TOP"))
dem_map  

# tmap::tmap_save(clusters, 'F:/Sync/PhD_Writing/Paper2/QGIS_map/Struct_Groups/cluster.svg')
# tmap::tmap_save(dem_map, 'F:/Sync/PhD_Writing/Paper2/QGIS_map/Struct_Groups/struct2.svg')
