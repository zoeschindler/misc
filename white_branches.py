########################################################################################################################
########################################################################################################################
# SOLVING THE WHITE BRANCH ISSUE
# ATTEMPT 1: RESAMPLE IMAGES BEFORE COLORING
########################################################################################################################
########################################################################################################################

# NOTES

# will probably have issues with flowers (overwriting flower pixels)
# good solution: depth images or flower classification to distinguish sky & flowers
# bad solution:  take as input whether the sky is blue or white

# opencv-python needs numpy to be installed
# opencv-python can be installed via pip

# TODO: function which reverses the renaming of the original images

########################################################################################################################

# PREPARATION

# import packages
import os
import glob
import cv2 as cv
import numpy as np
from random import sample


########################################################################################################################

# HELPER FUNCTIONS

# a) extract mean color value of pixels with 1
def window_mean(arr, threshold=1):
    if sum(arr[:, :, 3].flatten()) >= threshold:  # if there are at lest [threshold] non-sky pixels
        arr_new = arr[arr[:, :, 3] == 1][:, 0:3]  # get values of the non-sky pixels
        new_val = np.round(np.mean(arr_new, 0)).astype(int)  # get mean values per band
        return new_val
    else:
        return None


# b) extract color value of random pixel with 1
def window_random(arr, threshold=1):
    if sum(arr[:, :, 3].flatten()) >= threshold:  # if there are at lest [threshold] non-sky pixels
        arr_new = arr[arr[:, :, 3] == 1][:, 0:3]  # get values of the non-sky pixels
        new_val = arr_new[sample(range(arr_new.shape[0]), 1),]  # get random non-sky pixel
        return new_val
    else:
        return None


# c) extract color value of darkest (lowest L) pixel with 1
def window_darkest(arr, threshold=1):
    if sum(arr[:, :, 3].flatten()) >= threshold:  # if there are at lest [threshold] non-sky pixels
        arr_new = arr[arr[:, :, 3] == 1][:, 0:3]  # get values of the non-sky pixels
        new_val = arr_new[arr_new[:, 1] == min(arr_new[:, 1]),][0,]  # get the first pixel with the lowest lightness
        return new_val
    else:
        return None


########################################################################################################################

# MAIN FUNCTION

def branch_resample(img_path, old_path, out_path, kernel_size=3, threshold=1):
    # load image
    img_bgr = cv.imread(img_path)

    # convert to hsv
    # https://docs.opencv.org/4.x/df/d9d/tutorial_py_colorspaces.html
    img_hls = cv.cvtColor(img_bgr, cv.COLOR_BGR2HLS_FULL)

    # filter pixels, H: blue-ish, L: bright
    # values could be a bit fine-tuned
    lower_bound = np.array([130, 110, 0])  # ([130, 110, 100])
    upper_bound = np.array([170, 255, 255])
    mask = cv.inRange(img_hls, lower_bound, upper_bound)

    # bitwise conjunction with mask
    clipped = cv.bitwise_and(img_hls, img_hls, mask=mask)

    # show images
    # cv.imshow('img_hls', img_hls); cv.waitKey(0)
    # cv.imshow('mask', mask); cv.waitKey(0)
    # cv.imshow('clipped', clipped); cv.waitKey(0)

    # reverse mask values
    mask = mask - 255  # 1 = branch, 0 = sky

    # add mask as fourth channel
    img_hls_mask = np.dstack((img_hls, mask))

    # apply kernel_size x kernel_size convolution on the mask
    # https://docs.opencv.org/4.x/d4/d13/tutorial_py_filtering.html
    mask_blurred = cv.blur(mask.astype(float), (kernel_size, kernel_size))

    # get the position of all pixels: 0 < value < 1
    edges_indices = np.where((mask_blurred > 0) & (mask_blurred < 1))

    # loop through windows with the extracted locations as centers
    img_hls_new = img_hls.copy()
    buffer = int(kernel_size / 2)
    for i in range(len(edges_indices[0])):

        # get row & col of extracted location
        i_row = edges_indices[0][i]
        i_col = edges_indices[1][i]

        # get row & column boundaries
        row_lower = 0 if (i_row - buffer) < 0 else (i_row - buffer)
        row_upper = mask.shape[0] - 1 if (i_row + buffer) >= mask.shape[0] else (i_row + buffer)
        col_lower = 0 if (i_col - buffer) < 0 else (i_col - buffer)
        col_upper = mask.shape[1] - 1 if (i_col + buffer) >= mask.shape[1] else (i_col + buffer)

        # extract color & mask values
        window = img_hls_mask[row_lower:(row_upper + 1), col_lower:(col_upper + 1), ]

        # apply function to calculate new center pixel values
        # new_center = window_mean(window, threshold=1)
        # new_center = window_random(window, threshold=1)
        new_center = window_darkest(window, threshold=1)  # best

        # exchange center pixel values
        if new_center is not None:
            img_hls_new[i_row, i_col,] = new_center

    # reconvert to RGB color space
    img_bgr_new = cv.cvtColor(img_hls_new, cv.COLOR_HLS2BGR_FULL)

    # rename old image
    os.rename(img_path, old_path)

    # save image
    cv.imwrite(out_path, img_bgr_new)

    # return resampled image
    return img_bgr_new


########################################################################################################################

# LOOPING FUNCTION
def riscan_branch_resample(riscan_path, old_suffix="_old", out_suffix="", kernel_size=3, threshold=1):
    # get all scanning images
    scans_path = riscan_path + "/SCANS"
    img_paths = glob.glob(scans_path + "/ScanPos*/SCANPOSIMAGES/*.jpg", recursive=True)

    # loop through images
    for img_path in img_paths:
        # set paths
        old_path = img_path[:-4] + old_suffix + ".jpg"
        out_path = img_path[:-4] + out_suffix + ".jpg"

        # execute main function for all images
        branch_resample(img_path, old_path, out_path, kernel_size, threshold)


########################################################################################################################

# EXECUTION

# set paths
img_path = "C:/Daten/Arbeit/Testbilder/blueten.jpg"
old_path = img_path[:-4] + "_old.jpg"  # new name of input image
out_path = img_path[:-4] + "_out.jpg"  # new name of output image

# run function
branch_resample(img_path, old_path, out_path, 5)

# set path
riscan_path = "C:/Daten/Arbeit/2020-11-19test.RiSCAN"

# run looped function
riscan_branch_resample(riscan_path, "_old", "", 5)

########################################################################################################################
# ATTEMPT 2: RECOLOR POINTS WITH HIGH PLANARITY / LINEARITY
########################################################################################################################


########################################################################################################################
# ATTEMPT 3: RECOLOR POINTS CLOSE TO QSM & WITH HIGH PLANARITY / LINEARITY
########################################################################################################################


########################################################################################################################
