################################################################################
################################################################################
# RECUT SINGLETREE FROM POINT CLOUD
################################################################################
################################################################################

# recuts existing singletree point clouds from the stand point cloud
# -> useful if stand point cloud was slightly altered (rgb / filtering / ...)

################################################################################

# load packages
library(lidR)
library(future)

# set paths
catalog_path <- "D:/Walnuss/2021-12-15_iso1cm5k_devref15_octree1cmavg.las"
catalog_tile_path <- "D:/Walnuss/tiles/{ID}"
singletree_old_path <- "D:/Walnuss/singletree"
singletree_new_path <- "D:/Walnuss/singletree_clean"

# set variables
voxel_res <- 0.01 # should be similar to the voxel size used for downsampling
voxel_buffer <- voxel_res/2
cores <- 15 # number of available CPU cores

# set options
retile <- FALSE 
delete_tile <- FALSE

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
  ctg <- readTLSLAScatalog(catalog_path, chunk_size = 10, chunk_buffer = 0)
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
    all_z_vals <- unique(voxels$Z,2)
    for (z_val in all_z_vals) {
      
      # create raster
      z_subset <- voxels[voxels$Z == z_val, ]
      z_subset <- as.data.frame(z_subset)[, c("X", "Y", "V1")]
      # new_raster <- rasterFromXYZ(z_subset)
      # crs(new_raster) <- crs(las)
      new_raster <- raster(x = extent(
        min(z_subset$X) - voxel_buffer,
        max(z_subset$X) + voxel_buffer,
        min(z_subset$Y) - voxel_buffer,
        max(z_subset$Y) + voxel_buffer),
        crs=crs(las), resolution = voxel_res)
      
      # combine empty raster with data
      data_points <- SpatialPointsDataFrame(
        coords = z_subset[, c("X", "Y")],
        data = data.frame("V1" = z_subset[, "V1"]),
        proj4string = crs(las))
      # new_raster <- rasterize(data_points, new_raster, field="V1", update=TRUE)
      new_raster[cellFromXY(new_raster, data_points@coords)] <- data_points$V1

      # add raster values to point cloud
      las_z <- filter_poi(las, Z > (z_val - voxel_buffer) & Z <= (z_val + voxel_buffer))
      las_z <- merge_spatial(las_z, new_raster, "V1")
      
      # remove points outside voxels
      las_z <- filter_poi(las_z, V1 > 0)
      temp_las <- rbind(temp_las, las_z)
    }
    
    # return relevant points
    return(temp_las)
  }, .options = list(automerge = TRUE))
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