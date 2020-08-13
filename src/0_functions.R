############################################################################## #
## FUNCTIONS ###################################################################
############################################################################## #

library(gdalUtils)
library(dplyr)
library(tibble)
library(stars)


## TILING ######################################################################

tile <- function(file, tilesize, ...){
  # Provides information on tiles used to break down a raster.
  # file: path to raster file.
  # tilezie: maximum size of tiles. Either a numeric vector of length 1,
  # indicating both maximum width and height of the tiles, or a a numeric vector
  # of length 2, indicating first x-width and then y-height.

  file_info <- gdalinfo(file, raw_output = FALSE)

  tile_par <- list(x_total = file_info$columns,
                   y_total = file_info$rows,
                   tilesize_x = ifelse(length(tilesize) == 1,
                                       tilesize,
                                       tilesize[1]),
                   tilesize_y = ifelse(length(tilesize) == 1,
                                       tilesize,
                                       tilesize[2]))
  tile_par$n_x <- ceiling(tile_par$x_total / tile_par$tilesize_x)
  tile_par$n_y <- ceiling(tile_par$y_total / tile_par$tilesize_y)

  # Tidyverse methods are much faster here than apply. In theory, storing all
  # values as integers would make sense but R only supports 32-bit signed
  # integers -- maximum value given by .Machine$integer.max

  tiles <-
    cbind(tile_id = 1:(tile_par$n_x * tile_par$n_y),
          expand.grid(x_tile = (1:tile_par$n_x),
                      y_tile = (1:tile_par$n_y))) %>%
    as_tibble() %>%
    mutate(x_offset = (x_tile - 1) * tile_par$tilesize_x,
           y_offset = (y_tile - 1) * tile_par$tilesize_y,
           x_size = pmin((tile_par$x_total - x_offset),
                          tile_par$tilesize_x),
           y_size = pmin((tile_par$y_total - y_offset),
                                    tile_par$tilesize_y),
           id_offset = y_offset * tile_par$x_total + x_offset)

  return(tiles)
}


get_cell_ids <- function(tiles, tile_id, ...){
  # Returns cell IDs for a tile based on a dataframe containing tile information
  # (as created with `tile()`).

  id <- tile_id

  # Compute raster extent to avoid reading file info again
  raster_extent <- tiles %>%
    filter(x_offset == max(x_offset),
           y_offset == max(y_offset)) %>%
    transmute(x = x_offset + x_size,
              y = y_offset + y_size) %>%
    c(., recursive = TRUE)

  window <- tiles %>%
    filter(tile_id == id) %>%
    mutate(x_total = (id_offset - x_offset) / y_offset) %>%
    c(., recursive = TRUE)

  cell_ids <- 
    rep((1:window["x_size"]) + window["id_offset"], 
        window["y_size"]) +
    rep((0:(window["y_size"] - 1)) * raster_extent["x"], 
        each = window["x_size"]) %>%
    unname()

  return(cell_ids)
}

read_tile_stars <- function(file, tiles, tile_id, rename_band = NULL, ...){
  # Use the `read_stars` function to read a specified tile from a raster.
  # file: path to raster file.
  # tiles: overview of tiles a returned by `tile()`.
  # tile_id: the ID of the tile to be read, must reference a valid ID in `tiles`.
  # rename_band: name of the abdn holding the values. If `NULL`, the `stars`
  # default is used, i.e. the file name.
  tile <-
    read_stars(file,
               RasterIO = list(nXOff = tiles$x_offset[tile_id] + 1,
                               nYOff = tiles$y_offset[tile_id] + 1,
                               nXSize = tiles$x_size[tile_id],
                               nYSize = tiles$y_size[tile_id],
                               bands = 1))

  if(!is.null(rename_band)) {
      tile <- setNames(tile, rename_band)
    }

  cell_ids <- get_cell_ids(tiles, tile_id)

  tile <- tile %>%
    mutate(cell_id = cell_ids,
           tile_id = tile_id)

  return(tile)
}

execute_by_tile <- function(file, tilesize, 
                            exfun, arguments, 
                            progress = NULL, 
                            ...){
  # Executes a function over a set of raster tiles.
  # file: path to raster file.
  # tilezie: maximum size of tiles. Either a numeric vector of length 1,
  # indicating both maximum width and height of the tiles, or a a numeric vector
  # of length 2, indicating first x-width and then y-height.
  # exfun: name function to be executed for each tile.
  # arguments: list of arguments to be passed to `exfun`.
  # progress: Interval at which to send progress messages. Set to `NULL`
  # (default) to suppress progress messages.

  stack <- arguments
  stack$file <- file
  tiles <- tile(file = file, tilesize = tilesize)
  stack$tiles <- tiles

  results_list <- foreach(current_tile_id = tiles$tile_id, ...) %dopar% {
    if(!is.null(progress) == TRUE && current_tile_id %% progress == 0){
      print(paste("Processing tile", current_tile_id, "of", nrow(tiles)))
    }
    stack$tile_id <- current_tile_id
    stack$tile <- current_tile <- do.call(what = read_tile_stars, args = stack)
    tile_result <- do.call(what = exfun, args = stack)
    return(tile_result)
  } # End foreach loop
  return(results_list)
}

## CALCULATIONS ################################################################

area_by_group <- function(tile, fields = NULL, na_rm = FALSE, ...){
  # Computes the total area by group.
  # tile: a `stars` object for which to compute the area.
  # fields: vector or list of grouping variables used to summarize the
  # calculated areas. 
  # na_rm: Whether an area for cells containing NA should be computed and
  # returned. Setting to `TRUE` can decrease computation time if many NAs are
  # present.

  # Determine if NAs are to be removed before or after area calculations.
  if(na_rm){
    na_ratio <- sum(is.na(current_tile$lossperiod)) / prod(dim(current_tile))
    if(na_ratio >= 0.75){
      na_rm_before <- TRUE
    } else {
      na_rm_before <- FALSE}
  } else {
    na_rm_before <- na_rm
  }

  # Read tile and compute area
  tile_tbl <-
    st_as_sf(tile) %>%
    select_at(fields) %>%
    {if(na_rm && na_rm_before) na.omit(.) else .} %>%
    mutate(area_m2 = st_area(.)) %>%
    st_drop_geometry()

  tile_table <- group_by_at(tile_tbl, fields)

  tile_area <- tile_table %>%
    summarise(area_m2 = sum(area_m2), .groups = "drop") %>%
    {if(na_rm && !na_rm_before) na.omit(.) else .}

  return(tile_area)
}

## MISCELLANEOUS ###############################################################

split_geom_poly <- function(x){
  # Disaggregate geometry collection in polygon features and geometry
  # collections.
  # x: an `sf` object.
  x_split <- x %>%
  filter(st_is(., "GEOMETRYCOLLECTION")) %>%
    st_cast() %>%
    filter(st_is(., c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))) %>%
    rbind(filter(x, !st_is(x, "GEOMETRYCOLLECTION")))
  return(x_split)
}

get_id_centroids <- function(file, tiles, tile_id, ...){
  # Extracts centroids with cell IDs from a specified tile
  # file: path to raster file.
  # tiles: overview of tiles a returned by `tile()`.
  # tile_id: the ID of the tile to be read, must reference a valid ID in `tiles`.

  centroids <-
    read_tile_stars(file, tiles, tile_id) %>%
    select(id) %>%
    st_as_sf(as_points = TRUE, merge = FALSE)

  return(centroids)
}
