################################################################################
################################################################################
# USING TREEQSM IN R
################################################################################
################################################################################

# load packages
library(R.matlab)
library(rgl)
library(viridis)

# set path
mat_path <- "C:/Daten/Arbeit/Test_TreeQSM/test_qsm.mat"

################################################################################
# READ TREEQSM MAT-FILE - FUNCTIONS
################################################################################

get_as_df <- function(target_list, pattern="", dim=1, col=c(), preview=FALSE) {
  # loop through list entries matching the pattern
  # assumes that all selected list entries are scalars of the same length
  # transforms each list entry to a data frame column
  # dim is only used if the first list entry has multiple dimensions
  # dim = 1 if the data is stored as one variable per column
  
  # get list names
  list_names <- names(target_list)
  list_names <- list_names[grepl(pattern, list_names)]
  
  # prepare empty data frame
  first_entry = target_list[[list_names[1]]]
  if (length(first_entry) %in% dim(first_entry)) {
    list_nrow <- length(first_entry)
  } else {
    list_nrow <- dim(first_entry)[dim]
  }
  new_df <- data.frame(matrix(NA, nrow=list_nrow, ncol=0))
  
  # add data column-wise to the data frame
  for (list_name in list_names) {
    current_var <- data.frame(matrix(target_list[[list_name]], nrow=list_nrow)) # get data
    colnames(current_var) <- list_name
    new_df <- cbind(new_df, current_var)
  }
  
  # rename columns
  if (!is.null(col) & length(col) == ncol(new_df)) {
    colnames(new_df) <- col
    print("columns were successfully renamed")
  }
  
  # show head of the dataframe
  if (preview) print(head(new_df))
  
  # return data frame
  return(new_df)
}

################################################################################

read_qsm <- function(mat_path, qsm_var) {
  
  # read in data
  data_mat <- readMat(mat_path)
  data_mat <- data_mat[[qsm_var]][,,1]
  
  #####
  
  # extract input parameters
  input_parameters <- get_as_df((data_mat$rundata[,,1])$inputs[,,1])
  
  #####
  
  # prepare cylinder data
  cylinder_names_old <- names(data_mat$cylinder[,,1])
  cylinder_names_new <- c()
  for (idx in 1:length(cylinder_names_old)) {
    cylinder_names_new <- c(cylinder_names_new, ifelse(
      cylinder_names_old[idx] %in% c("start", "axis"),
      list(paste0(cylinder_names_old[idx], c("_X","_Y","_Z"))),
      cylinder_names_old[idx])[[1]])
  }
  
  # extract cylinder data
  cylinder <- get_as_df(data_mat$cylinder[,,1], dim=1, col=cylinder_names_new)
  
  #####
  
  # extract branch data
  branch <- get_as_df(data_mat$branch[,,1])
  
  #####
  
  # prepare treedata
  tree_mat <- data_mat$treedata[,,1]
  tree_names <- names(tree_mat)
  idx_loc <- which(tree_names == "location")
  tree_overview_mat <- tree_mat[1:(idx_loc-1)]
  tree_other_mat <- tree_mat[(idx_loc):length(tree_names)]
  
  # treedata - overview
  treedata_overview <- get_as_df(tree_overview_mat)
  
  # treedata - location
  location <- data.frame(
    "X" = tree_mat[[idx_loc]][1,1],
    "Y" = tree_mat[[idx_loc]][1,2],
    "Z" = tree_mat[[idx_loc]][1,3])
  
  # treedata - StemTaper
  stemtaper <- get_as_df(tree_other_mat, "StemTaper",
                         dim=2, col=c("distance_m","diameter_m"))
  stem_cylinders_xyz <- cylinder[
    cylinder$BranchOrder == 0,  # get stem cylinders
    c("start_X", "start_Y","start_Z","axis_X","axis_Y","axis_Z")]  # get coordinates
  stemtaper <- cbind(stemtaper, rbind(stem_cylinders_xyz, 0)) # add to stemtaper data frame
  
  # treedata - BranchOrder
  branchorder <- get_as_df(tree_other_mat, "BranchOrd")
  branchorder <- cbind(data.frame("BranchOrder" = 1:nrow(branchorder)), branchorder)
  
  # treedata - classes - all
  cylinder_dia <- get_as_df(tree_other_mat, "CylDia")  # diameter classes, 1cm
  cylinder_hei <- get_as_df(tree_other_mat, "CylHei")  # height classes, 1m
  cylinder_zen <- get_as_df(tree_other_mat, "CylZen")  # zenith classes, 10°
  cylinder_azi <- get_as_df(tree_other_mat, "CylAzi")  # azimuth classes, 10°
  
  # treedata - classes - branches
  branch_dia <- get_as_df(tree_other_mat, "BranchDia")  # diameter classes, 1cm
  branch_hei <- get_as_df(tree_other_mat, "BranchHei")  # height classes, 1m
  branch_zen <- get_as_df(tree_other_mat, "BranchZen")  # zenith classes, 10°
  branch_azi <- get_as_df(tree_other_mat, "BranchAzi")  # azimuth classes, 10°
  
  #####
  
  # prepare pmdistance data
  pmdist_mat <- data_mat$pmdistance[,,1]
  pmdist_names <- names(pmdist_mat)
  
  # extract pmdistance - overview
  pmdist_overview <- get_as_df(pmdist_mat[pmdist_names != "CylDist"])
  
  # extract pmdistamce - distances
  pmdist_distance <- data.frame("CylDist" = pmdist_mat[["CylDist"]])
  
  #####
  
  # return results
  return(list(
    "input_parameters" = input_parameters,
    "cylinder" = cylinder,
    "branch" = branch,
    "treedata_overview" = treedata_overview,
    "location" = location,
    "stemtaper" = stemtaper,
    "branchorder" = branchorder,
    "cylinder_dia" = cylinder_dia,
    "cylinder_hei" = cylinder_hei,
    "cylinder_zen" = cylinder_zen,
    "cylinder_azi" = cylinder_azi,
    "branch_dia" = branch_dia,
    "branch_hei" = branch_hei,
    "branch_zen" = branch_zen,
    "branch_azi" = branch_azi,
    "pmdist_overview" = pmdist_overview,
    "pmdist_distance" = pmdist_distance))
}

################################################################################
# READ TREEQSM MAT-FILE - EXECUTION
################################################################################

qsm <- read_qsm(mat_path, "QSM")

################################################################################
# PLOTTING QSMS - FUNCTIONS
################################################################################

plot_qsm <- function(data, col_var="BranchOrder", palette=turbo, light_scene=FALSE,
         bg_color="grey20", ax_color="white", window=c(500,700)) {
  # plots QSM cylinders in an rgl device
  # col_var:      which variable to use for coloring (branch / BranchOrder)
  # palette:      color palette to use (viridis, turbo, magma, ...)
  # light scene:  should the scene be lit
  # bg_color:     background color
  # ax_color:     axes color
  # window:       initial window size
  
  # extract cylinders
  if (is(data, "list")) {
    cylinder <- data$cylinder
  } else if (is(data, "data.frame")) {
    cylinder <- data
  } else {
    warning("input must be a QSM (list) or cylinders (data frame)")
    stop()
  }
  
  # calculate end points of cylinders
  cylinder$end_X = cylinder$start_X + cylinder$axis_X * cylinder$length
  cylinder$end_Y = cylinder$start_Y + cylinder$axis_Y * cylinder$length
  cylinder$end_Z = cylinder$start_Z + cylinder$axis_Z * cylinder$length
  
  # create color ramp
  cyl_vals <- unique(cylinder[,col_var])
  col_n <- length(cyl_vals)
  col_vec <- palette(col_n)
  
  # assign the colors to the cylinders
  cylinder$color <- NA
  for (idx in 1:col_n) {
    cylinder$color[cylinder[,col_var] == cyl_vals[idx]] <- col_vec[idx]
  }
  
  # plot the cylinders
  open3d()
  par3d(windowRect = c(50,50,window[1]+50,window[2]+50))
  bg3d(bg_color)
  for (i in 1:nrow(cylinder)) {
    cyl <- cylinder3d(
      center = cbind(
        c(cylinder$start_X[i], cylinder$end_X[i]), 
        c(cylinder$start_Y[i], cylinder$end_Y[i]),
        c(cylinder$start_Z[i], cylinder$end_Z[i])),
      radius = cylinder$radius[i],
      closed = -2)
    shade3d(cyl, col=cylinder$color[i], lit=light_scene)
  }
  axes3d(edges = c("x", "y", "z"), col = ax_color)
}

################################################################################
# PLOTTING QSMS - EXECUTION
################################################################################

plot_qsm(qsm, col_var="branch")
# or: plot_qsm(qsm$cylinder)

################################################################################