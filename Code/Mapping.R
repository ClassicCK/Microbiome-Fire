library(raster)
library(rgdal)
library(ggplot2)
library(cowplot)
library(broom)
library(ggmap)
library(sf)
library(basemaps)
library(gridExtra)

# Load in Fire Perimeter from Great Smoky Mountains National Park (previouly downloaded and provided in GitHub repo)
perimeter <- readOGR( 
  dsn = paste0(getwd(),"/data/GRSM_CT_FIRE_PERIMETER/"), 
  layer = "GRSM_CT_FIRE_FIRE_POLYGON",
  verbose = FALSE
)

st_crs(perimeter) #Check coordinate reference system
proj4string(perimeter) <- CRS("+init=epsg:26917") #Project to new CRS
perimeter <- spTransform(perimeter, CRS("+proj=longlat +datum=NAD83")) #Set CRS

fire <- st_as_sf(perimeter) #Create seperate shapefile set

# Load in Plots from NEON repository (previouly downloaded and provided in GitHub repo)
plots <- readOGR(
  dsn = paste0(getwd(),"/data/GRSM_NEON_PLOTS/"),
  layer = "IANDM_VEG_DBO_GRSM_NEON_PLOTS",
  verbose = FALSE
)

#Load in Burn Severity (previouly downloaded and provided in GitHub repo)
#burn <- read.csv(paste0(getwd(),"/data/GRSM_CT_FIRE_SOIL_BURN_SEVERITY.csv"))
severity <- readOGR( 
  dsn = paste0(getwd(),"/data/GRSM_CT_FIRE_SOIL_BURN_SEVERITY/"), 
  layer = "GRSM_CHIMNEY_TOPS_FIRE_SOIL_BURN_SEVERITY",
  verbose = FALSE
)

st_crs(severity) #Check coordinate reference system
proj4string(severity) <- CRS("+init=epsg:26917") #Project to new CRS
severity <- spTransform(severity, CRS("+proj=longlat +datum=NAD83 +init=epsg:26917")) #Set CRS

burn.severity <- st_as_sf(severity) #Create seperate shapefile set

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
register_google(key = "Place Key Here") #Will need Google API to run
tn_map <- get_map(location = c(lon = -83.475, lat = 35.70), zoom = 11, maptype = "roadmap") #Set bounding frame to state of Tennessee

p <- ggmap(tn_map) +
  geom_sf(data = fire, fill = "red", color="red", inherit.aes = FALSE, alpha = 0.2) +
  geom_sf(data = plots, color = c("#002244", "#002244", "#FB4F14", "#FB4F14",
                                 "#002244", "#FB4F14", "#FB4F14", 
                                 "#FB4F14", "#FB4F14"), inherit.aes = FALSE) +
  theme_cowplot()

#Save Plot
ggsave(filename= paste0(getwd(), "/Figures/Figure1A.pdf"), 
       plot = p, width = 180, height = 180, units=c("mm"), dpi=600) 

# Supplemental Figure with Burn Severity
p <- ggmap(tn_map) +
  geom_sf(data = burn.severity, aes(fill = burn.severity$Severity), inherit.aes = FALSE) +
  geom_sf(data = plots, color = c("#002244", "#002244", "#FB4F14", "#FB4F14",
                                  "#002244", "#FB4F14", "#FB4F14", 
                                  "#FB4F14", "#FB4F14"), inherit.aes = FALSE) +
  theme_cowplot()

#Save Plot
ggsave(filename= paste0(getwd(), "/Figures/FigureS5.pdf"), 
       plot = p, width = 180, height = 180, units=c("mm"), dpi=600) 
