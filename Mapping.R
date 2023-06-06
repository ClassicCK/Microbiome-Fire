library(raster)
library(rgdal)
library(ggplot2)

# Load in Fire Perimeter
perimeter <- readOGR( 
  dsn = paste0(getwd(),"/data/GRSM_CT_FIRE_PERIMETER/"), 
  layer = "GRSM_CT_FIRE_FIRE_POLYGON",
  verbose = FALSE
)

# Load in Plots
plots <- readOGR(
  dsn = paste0(getwd(),"/data/GRSM_NEON_PLOTS/"),
  layer = "IANDM_VEG_DBO_GRSM_NEON_PLOTS",
  verbose = FALSE
)


#Load in Burn Severity
burn <- read.csv(paste0(getwd(),"/data/GRSM_CT_FIRE_SOIL_BURN_SEVERITY.csv"))

#Subset to plots of interest
keep <- c("")