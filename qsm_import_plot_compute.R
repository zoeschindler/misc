################################################################################
################################################################################
# USING TREEQSM IN R
################################################################################
################################################################################

# load packages
library(R.matlab)
library(rgl)
library(viridisLite)
library(lidR)

# execute?
execute_example_1 <- FALSE
execute_example_2 <- TRUE

# CONTENT
# - READ TREEQSM MAT-FILE
# - PLOTTING QSM
# - CALCULATE QSM VIA MATLAB SERVER
# - EXAMPLE EXECUTION

################################################################################
# READ TREEQSM MAT-FILE
################################################################################

# loops through list entries matching a pattern
# assumes that all selected list entries are scalars of the same length
# transforms each list entry to a data frame column
get_as_df <- function(target_list, pattern="", dim=1, col=c(), preview=FALSE) {
  # target_list: list to be converted into a data frame
  # pattern:     pattern of the list entry names to be considered
  # dim:         only used if first entry of list has multiple dimensions,
  #              dim = 1 -> stored as one variable per column
  #              dim = 2 -> stored as one variable per row
  # col:         column names to override the existing names
  # preview:     print the first rows of the data frame?
  
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
  
  # add data column-wise to data frame
  for (list_name in list_names) {
    current_var <- data.frame(matrix(target_list[[list_name]], nrow=list_nrow)) # get data
    colnames(current_var) <- list_name
    new_df <- cbind(new_df, current_var)
  }
  
  # rename columns
  if (!is.null(col) & length(col) == ncol(new_df)) {
    colnames(new_df) <- col
  }
  
  # show head of data frame
  if (preview) print(head(new_df))
  
  # return data frame
  return(new_df)
}

################################################################################

# converts data from a matlab qsm to a list of dataframes
read_qsm <- function(data_in, qsm_var="QSM") {
  # data_in: path to a matlab file containing the qsm or the read in matlab file
  # qsm_var: name of the qsm in the matlab file
  
  # read in data
  if (is(data_in, "character")) {
    data_mat <- readMat(data_in)  # read matlab file from path
  } else if (is(data_in, "list")) {
    data_mat <- data_in  # already read matlab file
  } else {
    warning("input must be a path (character) or a qsm (list)")
    stop()
  }
  data_mat <- data_mat[[qsm_var]][,,1]
  
  # extract input parameters
  input_parameters <- get_as_df((data_mat$rundata[,,1])$inputs[,,1])
  
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
  
  # extract branch data
  branch <- get_as_df(data_mat$branch[,,1])
  
  # prepare treedata
  tree_mat <- data_mat$treedata[,,1]
  tree_names <- names(tree_mat)
  idx_loc <- which(tree_names == "location")
  tree_overview_mat <- tree_mat[1:(idx_loc-1)]
  tree_other_mat <- tree_mat[(idx_loc):length(tree_names)]
  
  # extract treedata - overview
  treedata_overview <- get_as_df(tree_overview_mat)
  
  #  extract treedata - location
  location <- data.frame(
    "X" = tree_mat[[idx_loc]][1,1],
    "Y" = tree_mat[[idx_loc]][1,2],
    "Z" = tree_mat[[idx_loc]][1,3])
  
  # extract treedata - StemTaper
  stemtaper <- get_as_df(tree_other_mat, "StemTaper",
                         dim=2, col=c("distance_m","diameter_m"))
  stem_cylinders_xyz <- cylinder[
    cylinder$BranchOrder == 0,  # get stem cylinders
    c("start_X", "start_Y","start_Z","axis_X","axis_Y","axis_Z")]  # get coordinates
  stemtaper <- cbind(stemtaper, rbind(stem_cylinders_xyz, 0)) # add to stemtaper data frame
  
  # extract treedata - BranchOrder
  branchorder <- get_as_df(tree_other_mat, "BranchOrd")
  branchorder <- cbind(data.frame("BranchOrder" = 1:nrow(branchorder)), branchorder)
  
  # extract treedata - classes - all
  cylinder_dia <- get_as_df(tree_other_mat, "CylDia")  # diameter classes, 1cm
  cylinder_hei <- get_as_df(tree_other_mat, "CylHei")  # height classes, 1m
  cylinder_zen <- get_as_df(tree_other_mat, "CylZen")  # zenith classes, 10°
  cylinder_azi <- get_as_df(tree_other_mat, "CylAzi")  # azimuth classes, 10°
  
  # extract treedata - classes - branches
  branch_dia <- get_as_df(tree_other_mat, "BranchDia")  # diameter classes, 1cm
  branch_hei <- get_as_df(tree_other_mat, "BranchHei")  # height classes, 1m
  branch_zen <- get_as_df(tree_other_mat, "BranchZen")  # zenith classes, 10°
  branch_azi <- get_as_df(tree_other_mat, "BranchAzi")  # azimuth classes, 10°
  
  # prepare pmdistance data
  pmdist_mat <- data_mat$pmdistance[,,1]
  pmdist_names <- names(pmdist_mat)
  
  # extract pmdistance - overview
  pmdist_overview <- get_as_df(pmdist_mat[pmdist_names != "CylDist"])
  
  # extract pmdistamce - distances
  pmdist_distance <- data.frame("CylDist" = pmdist_mat[["CylDist"]])
  
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
# PLOTTING QSM
################################################################################

# plots QSM cylinders in an rgl device
plot_qsm <- function(data, col_var="BranchOrder", palette=turbo, light_scene=FALSE,
         bg_color="grey20", ax_color="white", window=c(500,700)) {
  # col_var:     which variable to use for coloring (branch / BranchOrder)
  # palette:     color palette to use (viridis, turbo, magma, ...)
  # light scene: should the scene be lit
  # bg_color:    background color
  # ax_color:    axes color
  # window:      initial window size
  
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
  # (using apply might improve performance)
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
# CALCULATE QSM VIA MATLAB SERVER
################################################################################

# start a Matlab server
start_mat_server <- function(host="localhost", port=9999) {
  # resources on Matlab + R
  # https://mandymejia.com/2014/08/18/three-ways-to-use-matlab-from-r/
  # https://www.r-bloggers.com/2015/04/matlabr-a-package-to-calling-matlab-from-r-with-system/
  Matlab$startServer(port=port)
  server <- Matlab(host=host, port=port)
  isOpen <- open(server)
  if (!isOpen) R.oo::throw("MATLAB server is not running: waited 30 seconds.")
  return(server)
}

################################################################################

# close a Matlab server
stop_mat_server <- function(server) {
  close(server)
}

################################################################################

# read point cloud + convert coordinates to matrix
qsm_points <- function(path_points) {
  # path_points: path to a point cloud saved as las / laz / txt
  
  # get file extension
  file_split <- strsplit(path_points, split="[.]")[[1]]
  file_exten <- file_split[length(file_split)]
  
  # load data differently depending on file extension
  if (file_exten == "txt") {
    pts <- read.table(path_points)[,1:3]
    colnames(pts) <- c("X", "Y", "Z")
  } else if (file_exten %in% c("las","laz")) {
    pts <- readLAS(path_points)
    pts <- pts@data[,c("X","Y","Z")]
  } else {
    warning("read the coordinates manually and save them as matrix")
    stop()
  }
  
  # return point cloud as matrix with XYZ columns
  return(as.matrix(pts))
}

################################################################################

# get input values from Matlab server + change them
qsm_inputs <- function(server, changes) {
  # server:  running matlab server
  # changes: list containing the to be changed input parameters
  
  # get default input values
  evaluate(server,"clear inputs; create_input;")
  inputs <- getVariable(server, "inputs")
  inputs <- inputs$inputs[,,1]
  
  # change input values
  for (name in names(changes)) {
    if (name %in% names(inputs)) {
      inputs[[name]] <- changes[[name]]
    }
  }
  
  # convert back to matrices inside list
  for (name in names(inputs)) {
    inputs[[name]] <- matrix(inputs[[name]], nrow = 1)
  }
  
  # return mutated input values
  return(inputs)
}

################################################################################

# calculate QSM on Matlab server
qsm_treeqsm <- function(server, points, inputs, path_wd) {
  # server:  running matlab server
  # points:  matrix containing the point cloud in XYZ columns
  # inputs:  input parameter list (result of qsm_inputs)
  # path_wd: path where the necessary "results" folder will be created
  
  # create output path
  dir.create(file.path(path_wd, 'results'), showWarnings = FALSE)
  
  # set matlab working directory
  evaluate(server, paste0("cd ", path_wd, ";"))
  evaluate(server, "pwd")
  
  # turn off plots
  # does not always work though
  evaluate(server, "set(0,'DefaultFigureVisible','off');")
  evaluate(server, "figure('visible','off');")
  
  # set points in Matlab
  setVariable(server, points = points)
  evaluate(server,"points(1:5,:)")  # just to confirm the data is read in correctly
  
  # hand over input values to matlab
  setVariable(server, inputs = inputs)
  
  # run treeqsm
  evaluate(server, "QSM = treeqsm(points, inputs);")
  
  # pass QSM from MATLAB to R
  qsm <- getVariable(server, "QSM")
  qsm <- read_qsm(qsm, "QSM")
  
  # return qsm
  return(qsm)
}

################################################################################

# saving point cloud as mat file
# if the variable(s) should have specified names, use a list
qsm_points_to_mat <- function(server, point_matrix, path_mat) {
  # server:       running matlab server
  # point_matrix: points in a matrix with XYZ columns
  # path_mat:     path to the target matlab file, without extension
  
  # hand over variables to matlab
  setVariable(server, path_mat = path_mat)
  setVariable(server, point_matrix = point_matrix)
  
  # save point cloud in mat-file
  if (is(point_matrix, "list")) {
    # if point_matrix is a list of matrices
    evaluate(server, "point_matrix")
    evaluate(server, "save(path_mat, '-struct', 'point_matrix')")
  } else if (is(point_matrix, "matrix")) {
    # if point_matrix is a single matrix
    evaluate(server, "point_matrix(1:5,:)")  # just to confirm the data is read in correctly
    evaluate(server, "save(path_mat, 'point_matrix')")
  } else {
    warning("point_matrix must be a list of matrices or a single matrix")
    stop()
  }
}

################################################################################

# make_models_parallel / make_models
qsm_make_models <- function(server, path_wd, path_mat, inputs, qsm_name, qsm_num, parallel = T) {
  # server:   running matlab server
  # path_wd:  path where the necessary "results" folder will be created
  # path_mat: path to the target matlab file (result of qsm_points_to_mat)
  # inputs:   input parameter list (result of qsm_inputs)
  # qsm_name: name for the output file
  # qsm_num:  number of models per point cloud
  # parallel: parallel computing on / off
  
  # create output path
  dir.create(file.path(path_wd, 'results'), showWarnings = FALSE)
  
  # set matlab working directory
  evaluate(server, paste0("cd ", path_wd, ";"))
  evaluate(server, "pwd")
  
  # delete extension if it was 'accidentally' given
  path_mat <- strsplit(path_mat, split="[.]")[[1]][1]
  
  # hand over variables to matlab
  setVariable(server, path_mat = path_mat)
  setVariable(server, qsm_name = qsm_name)
  setVariable(server, qsm_num = qsm_num)
  setVariable(server, inputs = inputs)
  
  # compute the models
  if (parallel) {
    evaluate(server, "QSMs = make_models_parallel(path_mat, qsm_name, qsm_num, inputs);")
  } else {
    evaluate(server, "QSMs = make_models(path_mat, qsm_name, qsm_num, inputs);")
  }
  
  # pass QSMs from matlab to R
  evaluate(server, "QSM_n = max(size(QSMs))")
  QSM_n <- getVariable(server, "QSM_n")
  QSM_n <- QSM_n[[1]][1,1]
  QSMs <- list()
  for (idx in 1:QSM_n) {
    evaluate(server, paste0("currentQSM = QSMs(", idx, ")"))
    current_QSM <- getVariable(server, "currentQSM")
    QSMs[[idx]] <- read_qsm(current_QSM, "currentQSM")
  }
  
  # return QSMs
  return(QSMs)
}

################################################################################

# select optimal input parameters & QSM
qsm_select_optimum <- function(server, path_wd, path_mat, qsm_name, qsm_measure="all_mean_dis") {
  # server:      running matlab server
  # path_wd:     path where the necessary "results" folder will be created
  # path_mat:    path to the target matlab file (result of qsm_make_models)
  # qsm_name:    name for the output file
  # qsm_measure: measure used for selection
  
  # create output path
  dir.create(file.path(path_wd, "results"), showWarnings = FALSE)
  
  # set matlab working directory
  evaluate(server, paste0("cd ", path_wd, ";"))
  evaluate(server, "pwd")
  
  # hand over variables to matlab
  setVariable(server, path_mat = path_mat)
  setVariable(server, qsm_name = qsm_name)
  setVariable(server, qsm_measure = qsm_measure)
  
  # load previously run models
  # assumes qsms were saved with make_models_parallel (-> one variable stored in file)
  evaluate(server, "QSM_data = load(path_mat)")
  evaluate(server, "QSM_data_names = fieldnames(QSM_data)")
  evaluate(server, "QSMs = QSM_data.(QSM_data_names{1,1})")
  
  # select the best models
  evaluate(server, "[TreeData, OptModels, OptInputs, OptQSM] = select_optimum(QSMs, qsm_measure, qsm_name)")
  
  # pass results from matlab to R
  TreeData  <- getVariable(server, "TreeData")
  OptModels <- getVariable(server, "OptModels")
  OptInputs <- getVariable(server, "OptInputs")
  OptQSM    <- getVariable(server, "OptQSM")
  
  # reshape Treedata
  # out: list with data frames with attributes
  treedata_mat <- TreeData$TreeData[,,1]
  treedata_num <- ifelse(is.null(ncol(treedata_mat)), 1, ncol(treedata_mat))
  treedata_names <- unlist(ifelse(treedata_num == 1,
                                  treedata_mat$name,
                                  list(treedata_mat["name",1:treedata_num])))
  # if there is only one tree, add a dummy dimension
  if (treedata_num == 1) {
    treedata_mat <- as.list(cbind(treedata_mat, treedata_mat))
  }
  treedata_list <- list()
  # loop through columns
  for (col_idx in 1:treedata_num) {
    idx_loc <- which(names(treedata_mat[,col_idx]) == "location")
    col_name <- as.character(treedata_mat["name",col_idx])
    treedata_overview <- get_as_df(treedata_mat[1:(idx_loc-1),col_idx])
    treedata_overview$name <- col_name
    treedata_overview$value <- c("mean", "sd")
    treedata_list[[col_name]] <- treedata_overview
  }
  
  # reshape OptModels
  # out: list with lists with optimal models
  optmodels_mat <- OptModels$OptModels
  optmodels_list <- list()
  for (idx in 1:treedata_num) {
    best_input <- unlist(optmodels_mat[[idx]])
    best_QSM   <- unlist(optmodels_mat[[idx + treedata_num]])
    optmodels_list[[treedata_names[idx]]] <- list(
      "using_best_input" = best_input, "best_single_qsm"  = best_QSM)
  }
  
  # reshape OptInputs
  # out: list with data frames with input value 
  optinputs_mat <- OptInputs$OptInputs
  optinputs_values <- unlist(optinputs_mat)
  optinputs_length <- length(optinputs_values)/treedata_num
  optinputs_names <- unlist(dimnames(optinputs_mat))
  optinputs_list <- list()
  for (idx in 1:treedata_num) {
    optinputs_idx <- optinputs_values[(1:optinputs_length)+(idx-1)*optinputs_length]
    optinputs_idx <- data.frame(matrix(optinputs_idx, nrow = 1))
    colnames(optinputs_idx) <- optinputs_names
    optinputs_list[[treedata_names[idx]]] <- optinputs_idx
  }
  
  # reshape OptQSMs
  # out: list with all read in QSMs
  optqsm_mat <- OptQSM
  optqsm_list <- list()
  for (idx in 1:treedata_num) {
    evaluate(server, paste0("currentQSM = OptQSM(", idx, ")"))
    optqsm_idx <- getVariable(server, "currentQSM")
    optqsm_list[[treedata_names[idx]]] <- read_qsm(optqsm_idx, "currentQSM")
  }
  
  # return all reshaped outputs
  return(list(TreeData  = treedata_list,
              OptModels = optmodels_list,
              OptInputs = optinputs_list,
              OptQSMs   = optqsm_list))
}

################################################################################
# EXECUTION EXAMPLE 1
################################################################################

# set path
wd_path <- "D:/Test_TreeQSM"

# execute example
if (execute_example_1) {
  
  # start server
  mat_server <- start_mat_server()
  
  # read points
  point_matrix <- qsm_points(paste0(wd_path, "/tree.txt"))
  
  # set inputs
  input_list <- qsm_inputs(mat_server, list(name = "banane",
                                            Tria = 0,
                                            PatchDiam1 = 0.1,
                                            PatchDiam2Min = 0.02,
                                            PatchDiam2Max = 0.06,
                                            BallRad1 = 0.1 + 0.02,
                                            BallRad2 = 0.06 + 0.01,
                                            nmin1 = 5,
                                            savemat = 1,
                                            savetxt = 1))
  
  # calculate QSM
  qsm <- qsm_treeqsm(mat_server, point_matrix, input_list, wd_path)
  
  # close server
  close(mat_server)
  
  # read QSM from file
  qsm <- read_qsm(file.path(wd_path, "results", "QSM_banane_t1_m1.mat"), "QSM")
  
  # plot QSM
  plot_qsm(qsm)
}

################################################################################
# EXECUTION EXAMPLE 2
################################################################################

# set path
path_wd <- "D:/Test_TreeQSM"

# execute example
if (execute_example_2) {
  
  # start server
  server <- start_mat_server()
  
  # read multiple point clouds
  points_1 <- qsm_points(file.path(path_wd, "tree.txt"))
  points_2 <- qsm_points(file.path(path_wd, "tree_2.txt"))
  point_matrix <- list("apple" = points_1, "banana" = points_2)
  
  # save points in single mat file
  path_mat <- file.path(path_wd, "tree_both.mat")
  qsm_points_to_mat(server, point_matrix, path_mat)
  
  # set input parameters
  inputs <- qsm_inputs(server, list(Tria = 0,
                                    PatchDiam1 = c(0.10, 0.15),
                                    PatchDiam2Min = 0.02,
                                    PatchDiam2Max = 0.05,
                                    nmin1 = 5,
                                    BallRad1 = c(0.10, 0.15) + 0.02,
                                    BallRad2 = 0.05 + 0.01))
  
  # calculate multiple models
  # several models per input parameter combination + per point cloud
  multiple_models <- qsm_make_models(server, path_wd, path_mat, inputs,
                                     qsm_name = "multiple_models_both",
                                     qsm_num  = 2,
                                     parallel = TRUE)
  
  # choose optimal input parameter combinations + QSM per point cloud
  path_mat <- file.path(path_wd, "results", "multiple_models_both.mat")
  optmized_models <- qsm_select_optimum(server, path_wd, path_mat,
                                        qsm_name = "optimization_both",
                                        qsm_measure = "all_mean_dis")
  
  # stop server
  stop_mat_server(server)
}

################################################################################
