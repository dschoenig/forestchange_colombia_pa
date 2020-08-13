############################################################################## #
## DATA PREPARATION ############################################################
############################################################################## #

## PACKAGES AND FUNCTIONS ######################################################

library(foreach)
library(doParallel)
library(sf)
library(gdalUtils)
library(units)
library(lwgeom)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(stringr)

# Set working directory to project root
setwd("/home/danielschoenig/projects/forestchange_colombia_pa/")
# Source helper functions
source("src/0_functions.R")


############################################################################## #
## YEAR OF FOREST LOSS CATEGORIZED #############################################
############################################################################## #

# Classify forest loss according to pre (2013-2015) and post (2016-2018) peace
# treaty.


## PARAMETERS ##################################################################

# Threshold for forest cover in 2000 (in %).
fc_threshold <- 50

# Parent directory of GFC files;
# different datasets are assumed to reside in the corresponding subdirectories,
# e.g. all the lossyear granules are assumed to be in the directory "lossyear"
gfc_dir <- "data/gfc/"

# GFC version (including prefixed year)
version <- "2018-v1.6"

# GFC granules to act on
granules <-
  data.frame(lon = c("090W", "080W", 
                     "090W", "080W", "070W", 
                     "080W", "070W"),
              lat = c("20N", "20N", 
                      "10N", "10N", "10N", 
                      "00N", "00N"))

# Vector file containing country boundaries
boundary_file <- "data/colombia/gadm36_COL.gpkg"
# Name of the layer in the boundary file
boundary_layer <- "gadm36_COL_0"

# Output directory for processed tiles and mosaic
out_dir <- "results/lossyear_classified/"

# Path to save shapes of protected areas
protected_areas_file <- "results/protected_areas/protected_areas.gpkg"

# Path to save shapes of buffer zones
buffer_zones_file <- "results/protected_areas/buffer_zones.gpkg"

############################################################################## #
## CLASSIFICIATION #############################################################
############################################################################## #

# Generate file names for input and output

granules <- granules %>%
  mutate(lossyear = str_c(gfc_dir, "lossyear/",
                               "Hansen_GFC-", version,
                          "_lossyear_",
                          lat, "_", lon, ".tif", sep = ""),
         datamask = str_c(gfc_dir, "datamask/",
                          "Hansen_GFC-", version,
                          "_datamask_",
                          lat, "_", lon, ".tif", sep = ""),
         treecover2000 = str_c(gfc_dir, "treecover2000/",
                               "Hansen_GFC-", version,
                               "_treecover2000_",
                               lat, "_", lon, ".tif", sep = ""),
         classified = str_c(out_dir,
                            "lossyear_classified_",
                            str_pad(fc_threshold,
                                    width = 2, side = "left", pad = "0"),
                            "_", lat, "_", lon, ".tif", sep = "")
         )


# Loop over granules (will take some minutes)

cl_gfc <- makeCluster(6)
registerDoParallel(cl_gfc, 6)

foreach(current_granule = 1:nrow(granules),
        .packages = c("gdalUtils", "stars")) %dopar% {

          print(paste("Processing granule", current_granule, 
                      "of", nrow(granules), "..."))

          # Filename for tmp files
          tmp_raster_classified <- tempfile("classified_", fileext = ".tif")
          tmp_boundary_clipped <- tempfile("clipped_", fileext = ".gpkg")

          # First: Use GDAL raster calculator to reclassify along four 
          # conditions simultaneously:
          # (1) Exclude non-terrestrial landmass (i.e. where B != 1);
          # (2) Exclude forest cover below threshold value (fc_threshold)
          # (3) Exclude location where forest loss occured before 2013;
          # (4) reclassify the year of forest loss to periods before (2013-2015)
          # and after peace treaty (2016-2018).
          system(paste("gdal_calc.py",
                       "-A", granules$lossyear[current_granule],
                       "-B", granules$datamask[current_granule],
                       "-C", granules$treecover2000[current_granule],
                       paste0(" --outfile=", tmp_raster_classified),
                       "--calc='(B!=1)*255 +",
                       # Exclude non-terrestrial
                       paste0("(C<", fc_threshold, ")*(B==1)*255 +"),
                       # Exclude loss before 2013
                       paste0("(C>=", fc_threshold, ")*(B==1)*(A>=1)*(A<=12)*255 +"),
                       # No loss is 0
                       paste0("(C>=", fc_threshold, ")*(B==1)*(A==0)*0 +"),
                       # Before peace treaty is 1
                       paste0("(C>=", fc_threshold, ")*(B==1)*(A>=13)*(A<=15)*1 +"),
                       # After peace treaty is 2
                       paste0("(C>=", fc_threshold, ")*(B==1)*(A>=16)*(A<=18)*2'"),
                       "--NoDataValue=255 --type='Byte' --co='COMPRESS=LZW'"),
                 intern = TRUE)

          # Second: Reproject vector layer and cut to granule extent. If this
          # step is performed first, dimensions of output rasters will not match
          # the datamask and treecover rasters. Alternatively, `crop_to_cutline
          # = FALSE` may be used, but then the raster will cover the entire
          # extent of the original tiles (i.e. beyound the bounding box of the
          # boundary shape)

          # Info on raster SRS and extent
          granule_srs <- st_crs(read_stars(granules$lossyear[current_granule], 
                                           proxy = TRUE))
          granule_info <- gdalinfo(granules$lossyear[current_granule],
                                   raw_output = FALSE)
          granule_bbox <- paste(granule_info$bbox[1,1], granule_info$bbox[2,1],
                                granule_info$bbox[1,2], granule_info$bbox[2,2])

          # Reproject boundary to raster SRS and clip to raster bounding box.
          # Depending on the shape of the boundary, this can considerably reduce
          # processing time; e.g. factor 2 or more

          # Using gdalUtils::ogr2ogr() insteaddoes not work because clipdst is
          # duplicated. Use system call to ogr2ogr command line utility instead.
          ## ogr2ogr(src_datasource_name = boundary_file,
          ##         dst_datasource_name = tmp_boundary_clipped,
          ##         layer = boundary_layer,
          ##         # Change double to single quotes, to avoid parsing errors
          ##         t_srs = str_replace_all(granule_srs, "\"", "\'"),
          ##         clipdst = granule_bbox,
          ##         verbose = TRUE
          ##         )
          system(paste("ogr2ogr",
                       paste0("-t_srs EPSG:", granule_srs$epsg),
                       "-clipdst", granule_bbox,
                       tmp_boundary_clipped,
                       boundary_file, boundary_layer,
                       sep = " "), intern = TRUE)

          # Verify output
          ## ogrinfo(tmp_boundary_clipped, boundary_layer, so = TRUE)

          # Mask raster using the clipped boundary
          # This is the most time consuming step, up to 6 minutes per granule
          gdalwarp(srcfile = tmp_raster_classified,
                   dstfile = granules$classified[current_granule],
                   cutline = tmp_boundary_clipped,
                   cl = boundary_layer,
                   crop_to_cutline = TRUE,
                   dstnodata = 255,
                   ot = "Byte",
                   multi = TRUE,
                   wo = c("NUM_THREADS=ALL_CPUS"),
                   co = c("COMPRESS=LZW"), verbose = TRUE)

          # Verify output
          ## gdalinfo(tmp_raster_masked)

          file.remove(tmp_boundary_clipped, tmp_raster_classified)
          file.exists(granules$classified[current_granule])

        } # End loop over granules

stopCluster(cl_gfc)

# Create mosaic from tiles
lossyear_classified_file <- paste0(out_dir, "lossyear_classified_",
                                   str_pad(fc_threshold, 2, "left", "0"),
                                   "_all.tif")
mosaic_rasters(gdalfile = granules$classified[file.exists(granules$classified)],
               dst_dataset = lossyear_classified_file,
               co = c("COMPRESS=LZW"))


############################################################################## #
## PROTECTED AREAS #############################################################
############################################################################## #

# Load RUNAP info
protected_areas <- read_sf("data/colombia/runap2/runap2Polygon.shp")

# Filter records
protected_areas <- protected_areas %>%
  filter(categoria %in% c("Parque Nacional Natural",
                          "Reserva Natural")) %>%
  # Exclude areas before study period
  filter(fecha_act < as.Date("2013-01-01"))

# Replace Serrania de los Churumbelos with old shape (prtocetd area was extended
# in 2018):
# Load old shape
serrania_ch_2018 <- read_sf("data/colombia/Parques_Nacionales_Naturales_de_Colombia-shp/26c2779d-7287-4ffd-a230-7aca288032a9202045-1-lckxrr.r37op.shp")
# Find features
rid2018 <- which(serrania_ch_2018$Nombre == "SERRANÍA DE CHIRIBIQUETE")
rid2020 <- which(protected_areas$nombre == "Serrania de Chiribiquete")
# Replace geometry
protected_areas[rid2020, "geometry"] <- 
  st_transform(serrania_ch_2018[rid2018, "geometry"], 4686)

# Prepare for analysis
protected_areas <- protected_areas %>%
  select(nombre, categoria) %>%
  arrange(nombre) %>%
  rowid_to_column("pa_id") %>%
  mutate(area_total_m2 = st_area(geometry)) %>%
  st_transform(st_crs(lossyear_classified)$epsg) %>%
  transmute(pa_id, nombre, categoria,
            area_total_km2 = set_units(area_total_m2, km^2))

# Export features
write_sf(protected_areas, protected_areas_file)


############################################################################## #
## BUFFER ZONES ################################################################
############################################################################## #

# Create raw (i.e. overlapping) buffers
buffer_raw <- read_sf(protected_areas_file) %>%
  select(nombre, categoria) %>%
  arrange(nombre) %>%
  rowid_to_column("pa_id") %>%
  select(pa_id) %>%
  # Transform to projected CRS, unit: m
  st_transform(3116) %>%
  st_buffer(10000) %>%
  st_difference(st_union(st_transform(protected_areas, 3116))) %>%
  st_set_precision(1) %>%
  st_make_valid()


## SPLITTING OVERLAPPING BUFFERS ###############################################

# Load classified forest loss
lossyear_classified <- read_stars(lossyear_classified_file, proxy = TRUE)

# Identify and extract overlapping regions in raw buffers
buffer_overlaps <- buffer_raw %>%
  split_geom_poly() %>%
  arrange(pa_id) %>%
  st_intersection() %>%
  filter(n.overlaps > 1) %>%
  st_transform(st_crs(lossyear_classified)$epsg)

# Transform to projected CRS for nearest feature calculation and simplify to
# less than half of the width of a raster cell in the GFC data set (10m
# tolerance).
protected_areas_3116 <- protected_areas %>%
  st_transform(3116) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 10)

# Narrow down potential PAs for each overlapping region
overlaps_to_pa <- buffer_overlaps$origins

# List to hold results
overlaps_assigned_list <- list()

# Set up parallel backend
cl_ov <- makeCluster(6, outfile = "log/buffer_overlaps.log")
registerDoParallel(cl_ov, 6)

for(i in 1:nrow(buffer_overlaps)){
  overlap <-
    lossyear_classified %>%
    st_crop(buffer_overlaps[i,]) %>%
    st_as_stars() %>%
    mutate(cell_id = seq.int(1, prod(dim(.)))) %>%
    select(cell_id) %>%
    st_crop(buffer_overlaps[i,]) %>%
    st_as_sf(na.rm = TRUE) %>%
    st_transform(3116)
  
  # Reduce PA table to only include relevant ones (this saves a lot of
  # computation time)
  protected_areas_ov <- protected_areas_3116 %>%
    slice(overlaps_to_pa[[i]])
  
  # Seperate overlap into chunks for parellel execution
  chunks <- tibble(from = as.integer(seq(1, nrow(overlap), 1000)),
                   to = as.integer(c(from[c(2:length(from))]-1,
                                     nrow(overlap))))
  
  overlap_nearest_pa <-
    foreach(chunk_id = 1:nrow(chunks),
            .packages = c("sf", "stars", "dplyr"),
            .combine = rbind) %dopar% {
              if(chunk_id %% 10 == 1){
                print(paste("Processing chunk", chunk_id, "of", nrow(chunks),
                            "for overlap", i, "of", nrow(buffer_overlaps)))
              }
              chunk_result <-
                overlap[chunks$from[chunk_id]:chunks$to[chunk_id],] %>%
                mutate(nearest_id =
                         st_nearest_feature(., protected_areas_ov)) %>%
                group_by(nearest_id) %>%
                summarise(.groups = "drop") %>%
                transmute(pa_id = protected_areas_ov$pa_id[nearest_id])
              
              return(chunk_result)
            } # End foreach loop over chunks
  
  overlaps_assigned_list[[i]] <- overlap_nearest_pa %>%
    group_by(pa_id) %>%
    summarise(.groups = "drop")
  
  rm(overlap, overlap_nearest_pa)
} # End for loop over overlaps

stopCluster(cl_ov)

overlaps_assigned <- do.call(rbind, overlaps_assigned_list) %>%
  group_by(pa_id) %>%
  summarise(.groups = "drop") %>%
  rename(geom = geometry) %>%
  st_set_precision(1) %>%
  st_make_valid()


## FINALIZE BUFFER ZONES #######################################################

# Combine assigned overlaps with non-overlapping parts of buffer zones. This
# creates small triangular "holes" (i.e. inner rings). These holes cover less
# then half of a raster cell and will therefore not impact calculation of forest
# loss. For visualization and other purposes, the polygons were subsequently
# cleaned by hand.

buffer_zones <- buffer_raw %>%
  split_geom_poly() %>%
  arrange(pa_id) %>%
  st_intersection() %>%
  filter(n.overlaps == 1) %>%
  st_set_precision(1) %>%
  st_make_valid() %>%
  split_geom_poly() %>%
  select(pa_id) %>%
  rbind(overlaps_assigned) %>%
  group_by(pa_id) %>%
  summarise(.groups = "drop") %>%
  left_join(st_drop_geometry(protected_areas), by = "pa_id") %>%
  select(-area_total_km2) %>%
  mutate(area_total_m2 = st_area(geom)) %>%
  ## st_transform(st_crs(lossyear_classified)$epsg) %>%
  transmute(pa_id, nombre, categoria,
            area_total_km2 = set_units(area_total_m2, km^2)) %>%
  st_set_crs(3116)

# Export features
write_sf(buffer_zones, buffer_zones_file)

