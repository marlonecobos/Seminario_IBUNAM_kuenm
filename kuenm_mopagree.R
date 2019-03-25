#' Agreement of extrapolative areas of MOP layers
#'
#' @description kuenm_mopagree calculates raster layers that represent the agreement of strict
#' extrapolative areas among two or more climate models of a emission scenario in a
#' given time period. Various emission scenarios and time periods can be processed.
#'
#' @param mop.dir (character) name of the folder in which final models are (e.g., the output
#' folder after using the \code{\link{kuenm_mod}}) function. It is important to have only the folders
#' containing the models in this directory. It can be only one folder or multiple folders containing models
#' for the same species, created with distinct parameter settings. If models were projected, and the
#' distinct types of extrapolation were used, the name of the folders contained in this directory should
#' include a pattern describing the type of extrapolation used (e.g., "EC" for extrapolation and
#' clamping in Maxent).
#' @param in.format (character) format of model raster files. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "GTiff".
#' @param out.format (character) format of layers to be written in \code{out.dir}. Options are "ascii", "GTiff",
#' and "EHdr" = bil. Default = "GTiff".
#' @param current (character) if exist, pattern to look for when defining which is the scenario of current
#' projection to be excluded from calculations. If not defined, no current projection is assumed.
#' @param time.periods (character or numeric) pattern to be searched when identifying MOP layers for
#' distinct time projections. If not defined it is assumed that one time period was considered.
#' @param emi.scenarios (character) pattern to be searched for identifying distinct emission
#' scenarios (e.g., RCP). If not defined it is asumed that only one emission scenario was used.
#' @param out.dir (character) name of the output directory to be created in which subdirectories
#' containing raster layers of model statistics will be written. Default = "MOP_agremment".
#'
#' @return Folders named like the folders inside \code{mop.dir} containing raster layers in format
#' \code{out.format} that represent agreement of strict strapolative areas for each emission scenario
#' in a each time period. Folders will be written inside \code{out.dir}.
#'
#' @details
#' Users must be specific when defining the patterns that the function will search for. This patterns
#' must be part of the mop layer names so the function can locate each file without problems.
#' This function uses this system of work to avoid high demands of RAM while perfomring these analyses.
#'
#' @examples
#' # MOP layers must be already created before using this function.
#' # https://github.com/marlonecobos/kuenm
#'
#' # Arguments
#' mop_dir <- "MOP_results"
#' format <- "GTiff"
#' curr <- "current"
#' time_periods <- 2050
#' emi_scenarios <- c("RCP4.5", "RCP8.5")
#' out_dir <- "MOP_agremment"
#'
#' kuenm_mopagree(mop.dir = mop_dir, in.format = format, out.format = format,
#'                current = curr, time.periods = time_periods,
#'                emi.scenarios = emi_scenarios, out.dir = out_dir)


kuenm_mopagree <- function(mop.dir, in.format = "GTiff", out.format = "GTiff", current,
                           time.periods, emi.scenarios, out.dir = "MOP_agremment") {
  
  # installing needed packages if required
  pcakages <- c("raster", "rgdal")
  req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
  if (length(req_packages) > 0) {
    install.packages(req_packages, dependencies = TRUE)
  }
  
  cat("Preparing data for starting analyses, please wait...\n")
  
  if (!dir.exists(mop.dir)) {
    stop(paste(mop.dir, "does not exist in the working directory, check folder name",
               "\nor its existence."))
  }
  if (length(list.dirs(mop.dir, recursive = FALSE)) == 0) {
    stop(paste(mop.dir, "does not contain any subdirectory named as sets of projection variables;",
               "\neach subdirectory inside", mop.dir, "must containg at least one mop raster layer."))
  }
  if (missing(current)) {
    cat("Argument current is not defined, no current projection will be assumed.\n")
  }
  if (missing(time.periods)) {
    cat("Argument time.periods is not defined, only one time projection will be assumed.\n")
  }
  if (missing(emi.scenarios)) {
    cat("Argument emi.scenarios is not defined, only one emission scenario will be assumed.\n")
  }
  
  if (in.format == "ascii") {
    format <- ".asc$"
  }
  if (in.format == "GTiff") {
    format <- ".tif$"
  }
  if (in.format == "EHdr") {
    format <- ".bil$"
  }
  if (out.format == "ascii") {
    format1 <- ".asc"
  }
  if (out.format == "GTiff") {
    format1 <- ".tif"
  }
  if (out.format == "EHdr") {
    format1 <- ".bil"
  }
  
  # Reading mop names
  nstas <- list.files(mop.dir, pattern = paste0("MOP.*", format),
                      full.names = TRUE, recursive = TRUE)
  mopn <- list.files(mop.dir, pattern = paste0("MOP.*", format), recursive = TRUE)
  mopin <- unique(gsub("%.*", "%", mopn))
  mopin <- unique(gsub("^S.*/", "", mopin))
  
  # Folder for all outputs
  dir.create(out.dir)
  
  # Folders for sets
  sets <- dir(mop.dir)
  set_dirs <- paste0(out.dir, "/", sets)
  
  for (i in 1:length(sets)) {
    # Separating by sets
    ecl <- paste0(".*/", sets[i], "/.*")
    ecla <- gregexpr(ecl, nstas)
    eclam <- regmatches(nstas, ecla)
    setses <- unlist(eclam)
    
    # Folders per each set
    dir.create(set_dirs[i])
    
    # Copying current if exists
    if (!missing(current)) {
      cu <- paste0(".*", current, ".*")
      cur <- gregexpr(cu, setses)
      curr <- regmatches(setses, cur)
      curre <- unlist(curr)
      
      to_cur <- paste0(out.dir, "/", sets[i], "/", mopin, "_",
                       gsub(paste0(".*", current), current, curre))
      file.copy(from = curre, to =  to_cur)
      
    }
    
    # Time periods
    if (missing(time.periods)) {
      time.periods <- ""
      timep <- 1
      septi <- ""
    }else {
      timep <- time.periods
      septi <- "_"
    }
    
    for (j in 1:length(time.periods)) {
      # Separating by times if exist
      tp <- paste0(".*", time.periods[j], ".*")
      tpe <- gregexpr(tp, setses)
      tper <- regmatches(setses, tpe)
      tperi <- unlist(tper)
      
      ## Separating by scenarios if exist
      if (missing(emi.scenarios)) {
        emi.scenarios <- ""
        sepem <- ""
      } else {
        sepem <- "_"
      }
      
      ## If exist and more than one, separate by emission scenarios
      for (k in 1:length(emi.scenarios)) {
        ### Separating by scenarios if exist
        es <- paste0(".*", emi.scenarios[k], ".*")
        esc <- gregexpr(es, tperi)
        esce <- regmatches(tperi, esc)
        escen <- unlist(esce)
        
        ### Claculations
        a <- raster::stack(escen) # stack
        b <- raster::values(a) # matrix
        a <- a[[1]] # raster layer from stack
        
        dims <- dim(b) # dimensions of matrix
        b <- c(b) # matrix to vector
        
        #### Reclassify
        b[!is.na(b)] <- ifelse(na.omit(b) == 0, 1, 0)
        b <- matrix(b, dims) # vector to matrix again
        
        #### New layer with model agreement
        a[] <- apply(b, 1, sum)
        
        ### Writing files
        mopnams <- paste(out.dir, sets[i], paste0(mopin, septi, time.periods[j],
                                             sepem, emi.scenarios[k], "_agreement", format1), sep = "/")
        
        raster::writeRaster(a, filename = mopnams, format = out.format)
        
        cat(paste("   ", k, "of", length(emi.scenarios), "emission scenarios\n"))
      }
      cat(paste(" ", j, "of", length(time.periods), "time periods\n"))
    }
    cat(paste(i, "of", length(sets), "sets\n"))
  }
  cat(paste("\nCheck your working directory!!!", getwd(), "\n", sep = "    "))
}
