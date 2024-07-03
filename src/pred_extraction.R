# ------------------------------------------------------------------------
# Title: pred_extraction.R
# Description: R code to implement the creation and extraction of the 
#              covariables in a grid of Aragon.
# Authors : Jorge Castillo-Mateo and David Jiménez
# Date:
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Set environment Load required data
temp <- read.delim("data/text_files/temperature_aragon.txt")
stations <- read.delim("data/text_files/stations_info.txt")
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# PACKAGES
library("sf")
library("sp")
library("elevatr")
library("tidyverse")
library("rnaturalearth")
library("rnaturalearthdata")
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# MAIN

df <- st_as_sf(stations, coords = c("lon", "lat"),
  crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
colnames(df)[2] <- "elevation"

# Save longitude and latitude
df$x <- stations$lon
df$y <- stations$lat

# Project to adequate reference system
df <- st_transform(df, 2062)
df <- df[, c(1, 3, 4:5, 2)]
head(df)


# Obtain borders of Aragon
aragon <- ne_states(country = "Spain", returnclass = "sf")
aragon <- aragon[aragon$name %in% c("Zaragoza", "Huesca", "Teruel"),]
aragon <- st_union(aragon)
aragon <- st_transform(aragon, 2062)

# Build regular grid 10 x 10 km and intersect with Aragon
grid_aragon <- st_make_grid(aragon, cellsize = 10000, what = "centers")
grid_aragon <- st_intersection(grid_aragon, aragon)

# OBTAIN ELEVATION (in m)
grid_aragon <- data.frame(st_coordinates(grid_aragon))
colnames(grid_aragon) <- c("x", "y")
grid_aragon <- get_elev_point(grid_aragon, prj = "EPSG:2062", 
  src = "aws", z = 12, override_size_check = TRUE)
grid_aragon$elev_units <- NULL


# Recover longitude and latitude for the grid
grid_aragon$x <- st_coordinates(st_transform(grid_aragon, 
                                    "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))[,1]
grid_aragon$y <- st_coordinates(st_transform(grid_aragon, 
                                    "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))[,2]
grid_aragon <- grid_aragon[, c(1, 3:4, 2)]
head(grid_aragon)

# OBTAIN DISTANCE TO THE COAST (in km)
# Load the polygon of neighboring countries, join them, and projection 2062
world_map <- ne_countries(scale = "large", returnclass = 'sf')
european_union <- c("Algeria", "Andorra", "France", "Gibraltar", "Morocco", "Portugal", "Spain")
european_union_map <- 
  world_map %>% 
  filter(name %in% european_union)
background <- st_transform(european_union_map, 2062)
bbox_europe <- st_bbox(c(xmin = -10, ymin = 20, xmax = 50, ymax = 80), crs = st_crs(european_union_map))
european_union_map_cropped <- st_crop(european_union_map, bbox_europe)
european_union_map_cropped <- st_union(european_union_map_cropped)
european_union_map_cropped <- st_transform(european_union_map_cropped, 2062)

# Distance (in km) from each grid cell to the boundary
grid_aragon$dist <- as.vector(st_distance(
  st_boundary(european_union_map_cropped), 
  grid_aragon)) / 1000

# Distance (in km) from each station to the boundary
df$dist <- as.vector(st_distance(
  st_boundary(european_union_map_cropped), 
  df)) / 1000

# ------------------------------------------------------------------------
