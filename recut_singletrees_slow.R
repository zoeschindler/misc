################################################################################
################################################################################
# RECUT SINGLETREE FROM POINT CLOUD (SLOW)
################################################################################
################################################################################

# recuts existing singletree point clouds from the stand point cloud
# -> useful if stand point cloud was slightly altered (rgb / filtering / ...)

################################################################################

# load packages
library(lidR)
library(future)

# set paths
catalog_path <- "D:/Walnuss/2022-01-26_iso10cm5k_devref15_octree10cmavg.las"
catalog_tile_path <- "D:/Walnuss/tiles/{ID}"
singletree_old_path <- "D:/Walnuss/singletree_dummy"
singletree_new_path <- "D:/Walnuss/singletree_clean"

# set variables
voxel_res <- 0.03 # should be similar to the voxel size used for downsampling
voxel_buffer <- voxel_res/2
cores <- 15 # number of available CPU cores

# set options
retile <- FALSE
delete_tile <- FALSE
xyz_neighbours <- TRUE

# set z offset to include voxels below & above
z_offset <- ifelse(xyz_neighbours, voxel_res, 0)

# create output folders
dir.create(singletree_new_path, showWarnings = FALSE)
dir.create(basename(catalog_tile_path), showWarnings = FALSE)

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
  voxels <- voxel_metrics(tree, length(X), res = voxel_res, all_voxels = FALSE)
  
  # load relevant catalog tiles
  ctg <- readTLSLAScatalog(catalog_path, chunk_size = 10, chunk_buffer = 0.5)
  ctg <- catalog_intersect(ctg, extent(tree))
  
  # loop through catalog tiles
  plan(multisession, workers = cores)
  within <- catalog_apply(ctg, function(cluster) {
    
    # load tile
    las <- readLAS(cluster)
    if (is.empty(las)) return(NULL)
    
    # voxel subset
    voxel_subset <- voxels[
      voxels$X >= las@bbox["x","min"] - voxel_buffer & voxels$X <= las@bbox["x","max"] + voxel_buffer &
        voxels$Y >= las@bbox["y","min"] - voxel_buffer & voxels$Y <= las@bbox["y","max"] + voxel_buffer,]
    
    # setup empty las file
    temp_las <- add_attribute(las, 0, "V1")
    temp_las <- filter_poi(temp_las, Z == 1000)  # this should be an impossible criterium
    
    # create a raster for each voxel layer
    all_z_vals <- unique(voxels$Z)
    for (z_val in all_z_vals) {
      
      # get voxels of the selected z layers
      # z_subset <- voxels[voxels$Z == z_val, ]
      z_subset <- voxels[voxels$Z >= z_val-z_offset & voxels$Z <= z_val+z_offset, ]
      z_subset <- as.data.frame(z_subset)  # [, c("X", "Y", "V1")]
      
      # create empty list for raster storage
      raster_list <- list()
      
      # loop through vertical voxel layers (below, z_val, above)
      for (z_curr in unique(z_subset$Z)) {
        z_subset_curr <- z_subset[z_subset$Z == z_curr,]
        
        # create raster
        z_raster <- raster(x = extent(
          min(z_subset$X) - voxel_buffer,
          max(z_subset$X) + voxel_buffer,
          min(z_subset$Y) - voxel_buffer,
          max(z_subset$Y) + voxel_buffer),
          crs=crs(las), resolution = voxel_res)
        
        # combine empty raster with data
        data_points <- SpatialPointsDataFrame(
          coords = z_subset_curr[, c("X", "Y")],
          data = data.frame("V1" = z_subset_curr[, "V1"]),
          proj4string = crs(las))
        z_raster[cellFromXY(z_raster, data_points@coords)] <- data_points$V1
        
        # replace empty cells with 0
        z_raster <- reclassify(z_raster, cbind(NA, 0))
        
        # filter raster to include points nearby filled rasters
        if (xyz_neighbours) {
          z_raster <- raster::focal(z_raster, w = matrix(1/9, nrow = 3, ncol = 3), fun=max, pad = TRUE, padValue = 0)
        }
        
        # append raster list
        raster_list <- append(raster_list, z_raster)
      }
      
      # convert raster list to raster brick
      raster_multi_z <- stack(raster_list)
      
      # sum up all raster layers of raster brick
      if (nlayers(raster_multi_z) > 1) {
        z_raster <- calc(raster_multi_z, sum)
      }
      else {
        z_raster <- raster_multi_z[[1]]
      }
      
      # add raster values to point cloud
      las_z <- filter_poi(las, Z > (z_val - voxel_buffer) & Z <= (z_val + voxel_buffer))
      las_z <- merge_spatial(las_z, z_raster, "V1")
      
      # remove points outside voxels
      las_z <- filter_poi(las_z, V1 > 0)
      temp_las <- rbind(temp_las, las_z)
      
      # remove buffer
      temp_las <- filter_poi(temp_las, buffer == 0)
    }
    
    # return relevant points
    return(temp_las)
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
