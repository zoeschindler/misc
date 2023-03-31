################################################################################
################################################################################
# RISCAN SCAN POSITIONS TO CSV
################################################################################
################################################################################

# load packages
library(terra)

# set paths
in_path <- "C:/Users/Schindler/Downloads/ScanPos.kmz"
out_path <- "C:/Users/Schindler/Downloads"

# set ouput epsg
out_crs <- "epsg:32632"

################################################################################

# unzip kml
doc_path <- unzip(in_path, files = "doc.kml", exdir = out_path)

# read kml
points <- vect(doc_path)

# transform crs
crs(points) <- "epsg:4326"
points <- project(points, out_crs)

# get coordinates
x <- geom(points, df = T)$x
y <- geom(points, df = T)$y

# get names
name <- as.numeric(gsub("ScanPos", "", points$Name))

# get height
z <- as.numeric(sapply(strsplit(gsub("&nbsp", "", points$Description), ";"), "[[", 5))

# combine data
out <- data.frame("x" = x, "y" = y, "z" = z, "name" = name)

# save as csv
write.csv(out, file.path(out_path, "ScanPos.csv"), row.names = FALSE, quote = FALSE)

# delete kml
unlink(file.path(out_path, "doc.kml"))

################################################################################