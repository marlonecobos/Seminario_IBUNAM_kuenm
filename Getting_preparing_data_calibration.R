################
# Obtaining and organizing data for kuenm example
# BEFORE MODEL CALIBRATION
# Marlon E. Cobos, manubio13@gmail.com 
################

#########
# See repository https://github.com/marlonecobos/MNE_Curso_pract_calibration

#########
# Occurrence data
library(rgbif)

# new directory
dir.create("R/model_calib")
setwd("R/model_calib")

# getting the data from GBIF
species <- name_lookup(query = "Cynomys mexicanus",
                       rank="species", return = "data") # information about the species

# checking which taxon key returns information
for (i in 1:length(species$key)) {
  cat("key", (1:length(species$key))[i], "=",
      occ_count(taxonKey = species$key[i], georeferenced = TRUE), "\n")
}

key <- species$key[6] # using taxon key that returned information

occ <- occ_search(taxonKey = key, return = "data") # getting data using the taxon key

colnames(occ)[1:10] # checking names of columns to select the ones of interest

# keeping only georeferenced records
occurrences <- occ[!is.na(occ$decimalLatitude) & !is.na(occ$decimalLongitude),
                   c("name", "decimalLongitude", "decimalLatitude")]

# saving the set of occurrences 
dir.create("data")
write.csv(occurrences, "data/cmex_all.csv", row.names = FALSE)

# excluding records with no coordinates
occurrences <- occurrences[!is.na(occurrences$decimalLongitude) | !is.na(occurrences$decimalLatitude), ]

# excluding duplicates
occurrences$code <-  paste(occurrences$name, occurrences$decimalLongitude, # concatenating columns of interest
                           occurrences$decimalLatitude, sep = "_")

occurrences <- occurrences[!duplicated(occurrences$code), 1:4] # erasing duplicates

# excluding records with (0, 0) coordinates
occurrences <- occurrences[occurrences$decimalLongitude != 0 & occurrences$decimalLatitude != 0, 1:3]

# saving the new set of occurrences 
write.csv(occurrences, "data/cmex_clean.csv", row.names = FALSE)


# thinning
library(spThin)
thin(occurrences, lat.col = "decimalLatitude", long.col = "decimalLongitude", spec.col = "name",
     thin.par = 10, reps = 10, write.files = TRUE, max.files = 1, out.dir = "data", out.base = "cmex",
     write.log.file = FALSE, verbose = TRUE)


# training and testing data splitting. randomly 75% for training and 25% for testing
dir.create("C_mexicanus")
occ_thinn <- read.csv("data/cmex_thin1.csv")
occ_thinn$name <- gsub(" ", "_", occ_thinn$name)

all <- unique(occ_thinn)

all$check <- paste(all[,2], all[,3], sep = "_")
train <- all[sample(nrow(all), round((length(all[,1])/4 *3))), ]
test <- all[!all[,4] %in% train[,4], ]

all$check <- NULL
train$check <- NULL
test$check <- NULL

write.csv(all, "C_mexicanus/cmex_joint.csv", row.names = FALSE)
write.csv(train, "C_mexicanus/cmex_train.csv", row.names = FALSE)
write.csv(test, "C_mexicanus/cmex_test.csv", row.names = FALSE)



#########
# Bioclimatic layers
library(raster)
wc2.5min <- getData(name = "worldclim", var = "bio", res = 2.5)



#########
# Calibration area, M, or buffered points
library(sp)
library(rgeos)
library(rgdal)

WGS84 <- wc2.5min@crs # geographic projection

occ_sp <- SpatialPointsDataFrame(coords = occ_thinn[, 2:3], data = occ_thinn,
                                 proj4string = WGS84)

# planar projection
centroid <- gCentroid(occ_sp, byid = FALSE) # centroid of coordinates

## projecting with latitud and longitud in reference to centroid of occurrence points
AEQD <- CRS(paste0("+proj=aeqd +lat_0=", centroid@coords[2], " +lon_0=", centroid@coords[1], 
                   " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")) # planar projection

occ_pr <- spTransform(occ_sp, AEQD) # projection

# buffer
dist <- 200000 # distance for buffer in meters

buff_area <- gBuffer(occ_pr, width = dist) # buffer of 200 km

# reprojection to projection of raster layers
buff_arearp <- spTransform(buff_area, WGS84)

# masking bioclimatic layers
wc_masked <- mask(crop(wc2.5min, buff_arearp), buff_arearp) # cropping and masking layers
names(wc_masked)
sel_var <- wc_masked[[c(-8, -9, -18, -19)]] # excluding variables with artifacts

nums <- gsub("bio", "", names(sel_var))

dir.create("masked_variables")
writeRaster(sel_var, filename = "masked_variables/bio.tif", format = "GTiff", 
            bylayer = TRUE, suffix = nums)



#########
# PCA of variables for models
source("https://raw.githubusercontent.com/marlonecobos/ENM_manuals/master/Variables_processing/kuenm_rpca.R")

# no variable projections to distinct scenarios needed

var_folder <- "masked_variables" # name of folder with variables to be combined in distinct sets
out_folder <- "PCA_variables" # name of folder that will contain the sets 
in_format <- "GTiff" # other options available are "GTiff" and "EHdr" = bil 
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

pcs <- kuenm_rpca(vars.folder = var_folder, in.format = in_format, out.format = out_format, project = FALSE, 
                  n.pcs = n_pcs, out.dir = out_folder, return.in = TRUE)




#########
# Organizing variables for calibration
## M in C_mexicanus folder
dir.create("C_mexicanus/M_variables")

for (i in 1:3) {
  infolder <- paste0("C_mexicanus/M_variables/Set_", i)
  dir.create(infolder)
  
  pcnums <- 1:(dim(pcs[[3]])[3] - i + 1)
  writeRaster(pcs[[3]][[pcnums]], filename = paste0(infolder, "/pc.asc"), format = "ascii",
              bylayer = TRUE, suffix = pcnums)
}
