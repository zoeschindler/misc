################################################################################
################################################################################
# RECUT SINGLETREE FROM POINT CLOUD (FAST)
################################################################################
################################################################################

# recuts existing singletree point clouds from the stand point cloud
# -> useful if stand point cloud was slightly altered (rgb / filtering / ...)

################################################################################

# load packages
library(lidR)
library(future)

# set paths
catalog_path <- "D:/Walnuss/recut/2022-01-26_iso10cm5k_devref15_octree10cmavg.las"
catalog_tile_path <- "D:/Walnuss/recut/tiles/{ID}"
singletree_old_path <- "D:/Walnuss/recut/alt"
singletree_new_path <- "D:/Walnuss/recut/neu"

# set variables
voxel_res <- 0.01 # should be similar to the voxel size used for downsampling
voxel_buffer <- voxel_res/2
cores <- 15 # number of available CPU cores

# set options
retile <- FALSE
delete_tile <- FALSE
neighbours <- 1

# neighbour offset
offset <- voxel_res * neighbours

# create output folders
dir.create(singletree_new_path, showWarnings = FALSE)
dir.create(dirname(catalog_tile_path), showWarnings = FALSE)

################################################################################

# retile catalog
if (retile) {
  ctg <- readTLSLAScatalog(catalog_path, chunk_size = 10, chunk_buffer = 0)
  opt_output_files(ctg) <- catalog_tile_path
  plan(multisession, workers = cores)
  catalog_retile(ctg)
  plan(multisession, workers = 1)
}
catalog_path <- dirname(catalog_tile_path)

################################################################################

# get paths of all singletree point clouds
singletree_old_files <- list.files(singletree_old_path, pattern = "[.]las", full.names = TRUE)

# loop through all singletree point clouds
for (tree_file in singletree_old_files) {
  
  # load singletree 
  tree <- readTLSLAS(tree_file)
  
  # create voxels
  voxels <- voxel_metrics(tree, as.numeric(length(X) > 0), res = voxel_res, all_voxels = FALSE)
  
  # load relevant catalog tiles
  ctg <- readTLSLAScatalog(catalog_path, chunk_size = 10, chunk_buffer = 0.5) # TODO: use buffer + remove later
  ctg <- catalog_intersect(ctg, ext(tree))
  
  # loop through catalog tiles
  plan(multisession, workers = cores)
  within <- catalog_apply(ctg, function(cluster) {
    
    # load tile
    las <- readLAS(cluster)
    if (is.empty(las)) return(NULL)
    
    # only regard neighbours if specified
    if (neighbours > 0) {
      
      # duplicate voxels depending on neighbours
      rep_num <- length(-neighbours:neighbours)**3
      voxels_idx <- rep(1:nrow(voxels), each = rep_num)
      voxels_rep <- voxels[voxels_idx,]
      
      # make shift data
      shift <- expand.grid(X = seq(from = -offset, to = offset, by = voxel_res),
                           Y = seq(from = -offset, to = offset, by = voxel_res),
                           Z = seq(from = -offset, to = offset, by = voxel_res))
      
      # add shift to voxel dataset
      voxels_rep$X <- voxels_rep$X + shift$X
      voxels_rep$Y <- voxels_rep$Y + shift$Y
      voxels_rep$Z <- voxels_rep$Z + shift$Z
      
      # get voxels with an entry
      voxels <- unique(voxels_rep)
    }
    
    # add voxel attribute to point cloud
    vox <- lidR:::group_grid_3d(las@data$X, las@data$Y, las@data$Z, c(voxel_res, voxel_res), c(0, 0, 0.5 * voxel_res))
    las <- las |>
      add_lasattribute(vox[[1]], "x_vox", "x_vox") |>
      add_lasattribute(vox[[2]], "y_vox", "y_vox") |>
      add_lasattribute(vox[[3]], "z_vox", "z_vox")
    
    # merge point cloud with voxels (only keeps points with according voxel)
    las@data <- merge(las@data, voxels,
                      by.x = c("x_vox", "y_vox", "z_vox"),
                      by.y = c("X", "Y", "Z"))
    
    # remove voxel attributes
    las <- las |>
      remove_lasattribute('x_vox') |>
      remove_lasattribute('y_vox') |>
      remove_lasattribute('z_vox')
    
    # remove buffer
    las <- filter_poi(las, buffer == 0)
    
    # return relevant points
    if (nrow(las@data) == 0) {
      return(NULL)
    } else {
      return(las)
    }
  }, .options = list(automerge = TRUE, need_buffer = TRUE))
  plan(multisession, workers = 1)
  
  # save tree las file
  writeLAS(within, file.path(singletree_new_path, basename(tree_file)))
}

################################################################################

# delete tile folder
if (delete_tile) {
  unlink(catalog_path, recursive = TRUE) 
}

################################################################################
