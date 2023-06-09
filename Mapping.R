library(raster)
library(rgdal)
library(ggplot2)
library(cowplot)
library(broom)
library(ggmap)
library(sf)
library(basemaps)
library(gridExtra)

# Load in Fire Perimeter
fire <- sf::st_read(paste0(getwd(),"/data/GRSM_CT_FIRE_PERIMETER/GRSM_CT_FIRE_FIRE_POLYGON.shp"), package = "sf")
perimeter <- readOGR( 
  dsn = paste0(getwd(),"/data/GRSM_CT_FIRE_PERIMETER/"), 
  layer = "GRSM_CT_FIRE_FIRE_POLYGON",
  verbose = FALSE
)
st_crs(perimeter)
proj4string(perimeter) <- CRS("+init=epsg:26917")
perimeter <- spTransform(perimeter, CRS("+proj=longlat +datum=NAD83"))

fire <- st_as_sf(perimeter)

# Load in Plots
plots <- readOGR(
  dsn = paste0(getwd(),"/data/GRSM_NEON_PLOTS/"),
  layer = "IANDM_VEG_DBO_GRSM_NEON_PLOTS",
  verbose = FALSE
)

#Load in Burn Severity
burn <- read.csv(paste0(getwd(),"/data/GRSM_CT_FIRE_SOIL_BURN_SEVERITY.csv"))

#Subset to plots of interest
keep <- c("1", "2", "3", "7", "16", "55", "58", "59", "60")
sites <- plots[plots$LOC_NAME %in% keep,]
plots <- st_as_sf(sites)
plots <- plots[-4,]
plots <- plots[-4,]

#Add in Burn/Unburned Information
plots$fire <- c("Unburned", "Unburned", "Burned", "Burned", "Unburned", 
                "Burned", "Burned", "Burned", "Burned")

#Map Fire and Sites
register_google(key = "Place Key Here")
tn_map <- get_map(location = c(lon = -83.475, lat = 35.70), zoom = 11, maptype = "roadmap")

p <- ggmap(tn_map) +
  geom_sf(data = fire, fill = "red", color="red", inherit.aes = FALSE, alpha = 0.2) +
  geom_sf(data = plots, color = c("#002244", "#002244", "#FB4F14", "#FB4F14",
                                 "#002244", "#FB4F14", "#FB4F14", 
                                 "#FB4F14", "#FB4F14"), inherit.aes = FALSE) +
  theme_cowplot()

ggsave(filename= "/Users/christopher/Library/Mobile Documents/com~apple~CloudDocs/Documents/Duke/GRSM_Biome/Manuscript/Figure1A.pdf", 
       plot = p, width = 180, height = 180, units=c("mm"), dpi=600) 

# Try to create call out later
q <- ggplot() +
  borders(database = "state", reg = "tennessee") + 
  geom_rect(xmin = -83.7, xmax = -83.2, ymin = 35.5, ymax = 35.9, 
                                      fill = NA, colour = "black", size = 1.5) +
  theme_cowplot()
