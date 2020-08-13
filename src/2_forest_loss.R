############################################################################## #
## CALCULATE FOREST LOSS #######################################################
############################################################################## #

# Calculate the area of forest lost for protected areas, their buffer zones, and
# Colombia as a whole.

## PACKAGES AND FUNCTIONS ######################################################

library(foreach)
library(doParallel)
library(units)
library(lwgeom)
library(stars)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(stringr)


# Set working directory to project root
setwd("/home/danielschoenig/projects/forestchange_colombia_pa/")
# Source helper functions
source("src/0_functions.R")


## PARAMETERS ##################################################################


# Threshold for forest cover in 2000 (in %).
fc_threshold <- 50

# Output directory for results
forestchange_dir <- "results/forestchange/"

# Directory where to classified forest loss mosaics
lossyear_dir <- "results/lossyear_classified/"

# Path to protected areas file
protected_areas_file <- "results/protected_areas/protected_areas.gpkg"

# Path to buffer zones file
buffer_zones_file <- "results/protected_areas/buffer_zones.gpkg"


## LOAD FOREST LOSS DATA AND SHAPES ############################################

# Load forest loss data according to defined threshold
lossyear_classified_file <- paste0(lossyear_dir, "lossyear_classified_",
                                   str_pad(fc_threshold, 2, "left", "0"),
                                   "_all.tif")
lossyear_classified <- read_stars(lossyear_classified_file, proxy = TRUE)

# Load protected areas and buffer zones
protected_areas <- read_sf(protected_areas_file)
buffer_zones <- read_sf(buffer_zones_file)

############################################################################## #
## PROTECTED AREAS #############################################################
############################################################################## #



# Set up parallel backend
cl_area <- makeCluster(6, outfile = "log/forestchange_pa.log")
registerDoParallel(cl_area, 6)

# List to hold results
pa_summary_list <- list()

for(i in 1:nrow(protected_areas)){
  # Show progress
  print(paste("Processing protected area",
              i, "of", nrow(protected_areas)))

  # Clip raster and write to tmp file
  tmp_pa <- tempfile(pattern = "tmp_pa_", fileext = ".gpkg")
  tmp_loss_pa <- tempfile(pattern = "tmp_loss", fileext = ".tif")

  # Write shape to temporary file
  write_sf(protected_areas[i,], tmp_pa)

  # Clip raster
  gdalwarp(srcfile = lossyear_classified_file,
           dstfile = tmp_loss_pa,
           cutline = tmp_pa,
           crop_to_cutline = TRUE,
           dstnodata = 255,
           ot = "Byte",
           multi = TRUE,
           wo = c("NUM_THREADS=ALL_CPUS"),
           co = c("COMPRESS=LZW"),
           verbose = FALSE)
  # Pure R implementation is possible instead of GDAL, but is currently at least
  # an order of magnitude slower.

  # Calculate area for each period
  pa_result <-
    execute_by_tile(file = tmp_loss_pa,
                    tilesize = 250,
                    exfun = area_by_group,
                    arguments = list(rename_band = "lossperiod",
                                     fields  = c("tile_id", "lossperiod"),
                                     na_rm = TRUE),
                    progress = 100,
                    .export = c("get_cell_ids", "read_tile_stars"),
                    .packages = c("gdalUtils", "dplyr", "stars"))

  # Summarise and format output
  pa_summary <- pa_result %>%
    bind_rows %>%
    group_by(lossperiod) %>%
    summarise(area_m2 = sum(area_m2), .groups = "drop") %>%
    mutate(pa_id = protected_areas$pa_id[i],
           nombre = protected_areas$nombre[i],
           lossperiod = as.factor(lossperiod),
           area_m2 = set_units(area_m2, m^2)) %>%
    mutate(area_km2 = set_units(area_m2, km^2)) %>%
    select(pa_id, nombre, lossperiod, area_km2)
  pa_summary_list[[i]] <- pa_summary

  # Clean up tmp files
  file.remove(tmp_loss_pa, tmp_pa)
}

# Close cluster
stopCluster(cl_area)

# Summarise and prepare for export
forestchange_pa <-
  pa_summary_list %>%
  bind_rows() %>%
  mutate(lossperiod = factor(lossperiod)) %>%
  pivot_wider(names_from = lossperiod,
              values_from = area_km2,
              names_prefix = "lossperiod",
              values_fill = list(area_km2 = set_units(0, "km^2"))) %>%
  mutate(area_forest_km2 = lossperiod0 + lossperiod1 + lossperiod2) %>%
  rename(no_loss_km2 = lossperiod0,
         loss_1315_km2 = lossperiod1,
         loss_1618_km2 = lossperiod2) %>%
  mutate_at(vars(matches("km2")), set_units, quote(km^2))

# Export
forestchange_pa_file <- paste0(forestchange_dir, "forestchange_pa_",
                                   str_pad(fc_threshold, 2, "left", "0"),
                                   ".csv")
write_csv(forestchange_pa, forestchange_pa_file)


############################################################################## #
## BUFFER ZONES ################################################################
############################################################################## #


## CALCULATE FOREST LOSS #######################################################

# Set up parallel backend
cl_area <- makeCluster(6, outfile = "log/forestchange_bz.log")
registerDoParallel(cl_area, 6)

# List to hold results
bz_summary_list <- list()

for(i in 1:nrow(buffer_zones)){
  # Show progress
  print(paste("Processing buffer zone",
              i, "of", nrow(buffer_zones)))

  # Clip raster and write to tmp file
  tmp_bz <- tempfile(pattern = "tmp_bz_", fileext = ".gpkg")
  tmp_loss_bz <- tempfile(pattern = "tmp_loss", fileext = ".tif")


  # Write shape to temporary file
  write_sf(buffer_zones[i,], tmp_bz)

  # Clip raster
  gdalwarp(srcfile = lossyear_classified_file,
           dstfile = tmp_loss_bz,
           cutline = tmp_bz,
           crop_to_cutline = TRUE,
           dstnodata = 255,
           ot = "Byte",
           multi = TRUE,
           wo = c("NUM_THREADS=ALL_CPUS"),
           co = c("COMPRESS=LZW"),
           verbose = FALSE)
  # Pure R implementation is possible instead of GDAL, but is currently at least
  # an order of magnitude slower.

  # Calculate area for each period
  bz_result <-
    execute_by_tile(file = tmp_loss_bz,
                    tilesize = 250,
                    exfun = area_by_group,
                    arguments = list(rename_band = "lossperiod",
                                     fields  = c("tile_id", "lossperiod"),
                                     na_rm = TRUE),
                    progress = 100,
                    .export = c("get_cell_ids", "read_tile_stars"),
                    .packages = c("gdalUtils", "dplyr", "stars"))

  # Summarise and format output
  bz_summary <- bz_result %>%
    bind_rows() %>%
    group_by(lossperiod) %>%
    summarise(area_m2 = sum(area_m2), .groups = "drop") %>%
    mutate(pa_id = buffer_zones$pa_id[i],
           nombre = buffer_zones$nombre[i],
           lossperiod = as.factor(lossperiod),
           area_m2 = set_units(area_m2, m^2)) %>%
    mutate(area_km2 = set_units(area_m2, km^2)) %>%
    select(pa_id, nombre, lossperiod, area_km2)
  bz_summary_list[[i]] <- bz_summary

  # Clean up tmp files
  file.remove(tmp_loss_bz, tmp_bz)
}

# Close cluster
stopCluster(cl_area)

# Summarise and prepare for export
forestchange_bz <-
  bz_summary_list %>%
  bind_rows() %>%
  mutate(lossperiod = factor(lossperiod)) %>%
  pivot_wider(names_from = lossperiod,
              values_from = area_km2,
              names_prefix = "lossperiod",
              values_fill = list(area_km2 = set_units(0, "km^2"))) %>%
  mutate(area_forest_km2 = lossperiod0 + lossperiod1 + lossperiod2) %>%
  rename(no_loss_km2 = lossperiod0,
         loss_1315_km2 = lossperiod1,
         loss_1618_km2 = lossperiod2) %>%
  mutate_at(vars(matches("km2")), set_units, quote(km^2))

# Export
forestchange_bz_file <- paste0(forestchange_dir, "forestchange_bz_",
                               str_pad(fc_threshold, 2, "left", "0"),
                               ".csv")
write_csv(forestchange_bz, forestchange_bz_file)


############################################################################## #
## COLOMBIA ####################################################################
############################################################################## #

# Set up parallel backend
cl_area <- makeCluster(6, outfile = "log/forestchange_col.log")
registerDoParallel(cl_area, 6)

# Calculate forest loss
forestchange_col_summary <-
  execute_by_tile(file = lossyear_classified_file,
                  tilesize = 250, 
                  exfun = area_by_group,
                  arguments = list(rename_band = "lossperiod",
                                   fields  = c("tile_id", "lossperiod"),
                                   na_rm = TRUE),
                  progress = 100,
                  .export = c("get_cell_ids", "read_tile_stars"),
                  .packages = c("gdalUtils", "dplyr", "stars"))

# Close cluster
stopCluster(cl_area)

# Summarise and prepare for export
forestchange_col <-
  forestchange_col_summary %>%
  bind_rows() %>%
  group_by(lossperiod) %>%
  summarise(area_m2 = sum(area_m2), .groups = "drop") %>%
  transmute(lossperiod,
            area_km2 = set_units(area_m2, km^2)) %>%
  mutate(lossperiod = factor(lossperiod)) %>%
  pivot_wider(names_from = lossperiod, 
              values_from = area_km2, 
              names_prefix = "lossperiod",
              values_fill = list(area_km2 = set_units(0, "km^2"))) %>%
  mutate(area_forest_km2 = lossperiod0 + lossperiod1 + lossperiod2) %>%
  rename(no_loss_km2 = lossperiod0,
         loss_1315_km2 = lossperiod1,
         loss_1618_km2 = lossperiod2) %>%
  mutate_at(vars(matches("km2")), set_units, quote(km^2))

# Export
forestchange_col_file <- paste0(forestchange_dir, "forestchange_col_",
                               str_pad(fc_threshold, 2, "left", "0"),
                               ".csv")
write_csv(forestchange_col, forestchange_col_file)
