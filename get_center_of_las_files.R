################################################################################
################################################################################
# GET CENTERS OF MULTIPLE LAS FILES
################################################################################
################################################################################

# load packages
library(lidR)
library(terra)

# set paths
in_path <- "C:/Users/Schindler/Downloads/in"
out_path <- "C:/Users/Schindler/Downloads/out.shp"

# get all file paths
files <- list.files(in_path, pattern = "[.]las", full.names = TRUE)

# prepare storage
x_all <- c()
y_all <- c()
name_all <- c()

for (file in files) {
  
  # read las file
  las <- readLAS(file)
  
  # get name
  name_curr <- strsplit(basename(file), "[.]")[[1]][1]
  
  # get center x and y coordinate
  x_curr <- mean(c(ext(las)$xmin, ext(las)$xmax))
  y_curr <- mean(c(ext(las)$ymin, ext(las)$ymax))
  
  # save in list
  name_all <- c(name_all, name_curr)
  x_all <- c(x_all, x_curr)
  y_all <- c(y_all, y_curr)
}

# convert to vector
centers <- vect(matrix(c(x_all,y_all), ncol = 2), atts = data.frame("name" = name_all), crs = crs(las)@projargs)

# save vector as shp
writeVector(centers, out_path, overwrite = T)

################################################################################