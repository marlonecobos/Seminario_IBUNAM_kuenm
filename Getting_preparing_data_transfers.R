################
# Obtaining and organizing data for using kuenm
# AFTER MODEL CALIBRATION
# Marlon E. Cobos, manubio13@gmail.com 
################

#########
# Data and organization used for model calibration was obtained using 
# the script Getting_preparing_data_calibration.R in the repository
# https://github.com/marlonecobos/Seminario_IBUNAM_kuenm

#########
# Working directory
setwd("R/model_calib")

#########
# Bioclimatic layers for future scenarios
library(raster)

wc2.5min <- getData(name = "worldclim", var = "bio", res = 2.5)
rcp45cc <- getData(name = "CMIP5", var = "bio", res = 2.5, rcp = 45, model = "CC", year = 50)
rcp85cc <- getData(name = "CMIP5", var = "bio", res = 2.5, rcp = 85, model = "CC", year = 50)
rcp45mc <- getData(name = "CMIP5", var = "bio", res = 2.5, rcp = 45, model = "MC", year = 50)
rcp85mc <- getData(name = "CMIP5", var = "bio", res = 2.5, rcp = 85, model = "MC", year = 50)


# Masking to area of interest
library(maps)

## getting mexico shapefile from maps
WGS84 <- wc2.5min@crs # geographic projection

mexico_map <- maps::map(database = "world", region = "Mexico", fill = TRUE, plot = FALSE) # map of the world

mexico_pol <- sapply(strsplit(mexico_map$names, ":"), function(x) x[1]) # preparing data to create polygon
mexico <- maptools::map2SpatialPolygons(mexico_map, IDs = mexico_pol, proj4string = WGS84) # map to polygon

## masking bioclimatic layers to mexico
mx_masked <- lapply(list(wc2.5min, rcp45cc, rcp85cc, rcp45mc, rcp85mc), 
                    function(x){
                      x <- x[[c(-8, -9, -18, -19)]] # selecting layers
                      x <- crop(x, mexico) # cropping
                      mask(x, mexico) #  masking 
                    })

nums <- gsub("bio", "", names(mx_masked[[1]]))
scenarios <- c("current", "cc_45", "cc_85", "mc_45", "mc_85")

infolder <- "masked_transfer_variables"
dir.create(infolder)

lapply(1:5, function(x){
  dir.create(paste0(infolder, "/", scenarios[x]))
  writeRaster(mx_masked[[x]], filename = paste0(infolder, "/", scenarios[x], "/bio.tif"), format = "GTiff", 
              bylayer = TRUE, suffix = nums)
})

#########
# PCA of variables for transfering models
source("https://raw.githubusercontent.com/marlonecobos/ENM_manuals/master/Variables_processing/kuenm_rpca.R")

# no variable projections to distinct scenarios needed

var_folder <- "masked_variables" # name of folder with variables to be combined in distinct sets
proj_folder <- "masked_transfer_variables" # name of the folder containing one or more folders with variables for other scenarios
out_folder <- "PCA_variables_transfer" # name of folder that will contain the sets 
in_format <- "GTiff" # other options available are "GTiff" and "EHdr" = bil 
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

pcs <- kuenm_rpca(vars.folder = var_folder, in.format = in_format, out.format = out_format, project = TRUE, 
                  proj.vars = proj_folder, n.pcs = n_pcs, out.dir = out_folder, return.in = TRUE)

pcs <- pcs[paste0("PCRasters_", scenarios)] # in order
  
  
#########
# Organizing variables for transfers
## G in C_mexicanus folder
dir.create("C_mexicanus/G_variables")

selected_set <- as.character(read.csv("C_mexicanus/Calibration_results/best_candidate_models_OR_AICc.csv")[1, 1]) 
selected_set <- strsplit(selected_set, "_")[[1]][5:6]
selected_set <- paste(selected_set, collapse = "_")

pcs_selected_set <- list.files(path = "C_mexicanus/M_variables/Set_2/", pattern = ".asc$")
pcnums <- as.numeric(gsub(".asc", "", gsub("pc_", "", pcs_selected_set)))
  
Gset_folder <- paste0("C_mexicanus/G_variables/", selected_set)
dir.create(Gset_folder)


lapply(1:5, function(i){
  infolder <- paste0(Gset_folder, "/", scenarios[i])
  dir.create(infolder)
  
  writeRaster(pcs[[i]][[pcnums]], filename = paste0(infolder, "/pc.asc"), format = "ascii",
              bylayer = TRUE, suffix = pcnums)
})
