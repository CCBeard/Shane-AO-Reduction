import os
import numpy as np
import scipy
#import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from scipy.ndimage import interpolation as interp

from skimage.registration import phase_cross_correlation



## function to plot an image cube
## for this code a "cube" means a stack of image data arrays

def plot_grid(datacube,imagenames):
    no_A = len(datacube) ## number of individual images in the cube
    xplots = int(np.around(np.sqrt(no_A))) ## number of image grid columns
    yplots = xplots + 1 ## number of image grid rows--sometimes there are one fewer, but that's okay

#     print no_A, xplots, yplots ## this line is for troubleshooting

    gs = gridspec.GridSpec(yplots, xplots) ## define the image grid
    plt.figure(figsize=(15,15)) ## set the figure size
    for i in range(no_A):
        ## plots each individual image within the image grid:
        B = datacube[i]
        plt.subplot(gs[i])
        plt.imshow(np.log10(B), origin='lower', cmap='gray')
        plt.title(imagenames[i])


#############################
#Function for plotting

def plot_grid_background(datacube,imagenames):
    no_A = len(datacube) ## number of individual images in the cube
    xplots = int(np.around(np.sqrt(no_A))) ## number of image grid columns
    yplots = xplots + 1 ## number of image grid rows--sometimes there are one fewer, but that's okay

#     print no_A, xplots, yplots ## this line is for troubleshooting

    gs = gridspec.GridSpec(yplots, xplots) ## define the image grid
    plt.figure(figsize=(15,15)) ## set the figure size
    for i in range(no_A):
        ## plots each individual image within the image grid:
        B = datacube[i]
        plt.subplot(gs[i])
        plt.imshow(np.log10(B), origin='lower', cmap='gray',vmin=2,vmax=3)
        plt.title(imagenames[i])




#Function for identifying centroid position and counts
from scipy.ndimage.filters import median_filter, gaussian_filter
def guess_gaussian_parameters(d):
    """
    This image guesses the maximum intensity of the
    image by obtaining a smoothed version of the image
    via median + gaussian filtering. Then, finds the
    maximum of the smoothed image. Using the median-filtered
    image, the algorithm also estimates the width of the PSF and
    with this value and an estimate of the volume, the amplitude of
    such gaussian.

    Input

    d       Numpy array that contains the values of the pixels.

    Output

    x0      x coordinate of the image maximum

    y0      y coordinate of the image maximum

    sigma   Width of the PSF

    A       Estimated amplitude of the gaussian function
    """

    # First, smooth the image with a median filter. For this, find
    # the optimal window as the square-root of the geometric mean
    # between the sizes of the array. This step is good to kill outliers:
    window = int(np.sqrt(np.sqrt(d.shape[0]*d.shape[1])))
    if window % 2 == 0:
        window = window + 1


    # Originally there was a median filter here, but it really adds a lot of time to the routine
    #d_mfiltered = median_filter(d,size = window)
    # Going to try without it to see if it makes a significant difference for vetoing
    d_gfiltered = gaussian_filter(d,sigma = window)

    # Now, find the maximum of this image:
    y0, x0 = np.where(d_gfiltered == np.max(d_gfiltered))

    # Take the first element. This helps in two cases: (1) only one maximum has
    # been found, the outputs are numpy arrays and you want to extract the numbers
    # and (2) in case there are several maximums (this has not happened so
    # far but could in principle), pick the first of them:
    y0 = y0[0]
    x0 = x0[0]

    # Now estimate the width of the PSF by taking a "cross-plot" using the
    # maximum values found:
    x_cut = d[:,int(x0)]
    sigma_x = (np.sum(x_cut*(np.abs(np.arange(len(x_cut))-y0)))/np.sum(x_cut))/3.
    y_cut = d[int(y0),:]
    sigma_y = (np.sum(y_cut*(np.abs(np.arange(len(y_cut))-x0)))/np.sum(y_cut))/3.
    sigma = np.sqrt(sigma_x*sigma_y)

    # (Under) estimate amplitude assuming a gaussian function:
    A = (np.sum(d-np.median(d))/(2.*np.pi*sigma**2))

    return x0,y0,sigma,2.*A

#function for sigma clipping

def get_neighbors(i,j,data):
        #make an array of the values of all the neighbors of this coordinate

        #i: y coordinate that we are looking at
        #j: x coordinate we are looking at
        neighbors = []
        try:
            neighbors.append(data[i+1][j])
            neighbors.append(data[i-1][j])
            neighbors.append(data[i][j+1])
            neighbors.append(data[i][j-1])
            neighbors.append(data[i+1][j+1])
            neighbors.append(data[i-1][j-1])
            neighbors.append(data[i-1][j+1])
            neighbors.append(data[i+1][j-1])
            return neighbors

        except IndexError:
            return data[i][j]

def sigma_clip(object_list, flat_darkcor_data_out):

    flat_darkcor_sigmacut_data = {}
    skyframe_data = {}


    for img in object_list:
        data = flat_darkcor_data_out[img][450:1050,850:1450]
        print(img)
        #create a new array
        data2 = np.zeros((600,600))

        for i in range(len(data)):
             for j in range(len(data[0])):
                neighbors = get_neighbors(i,j,data)
                nvals = np.median(neighbors)
                if np.abs(nvals-data[i][j]) < 10:
                    data2[i][j] = data[i][j]
                else:
                    data2[i][j] = nvals

        flat_darkcor_sigmacut_data[img] = data2

    return flat_darkcor_sigmacut_data


def image_shift(object_list, center, sky_flat_darkcor_data_out, datadir):
    #Want to write something that checks all of these and finds the one with the highest S/N ratio
    #I'm thinking we can take the ratio of the centroid counts to some random background area
    #Will have to use that guess Gaussian parameters to find the centroid
    highest = 0
    zero_shift_image = center[0]
    for i in range(len(center)):
        data = sky_flat_darkcor_data_out[center[i]]

        x0,y0,sigma,A = guess_gaussian_parameters(data)
        centroid_window = data[y0-1:y0+1,x0-1:x0+1]
        cent_counts = np.mean(centroid_window)
        x1 = x0-200
        y1 = y0-200
        noise_window = data[y1-5:y1+5,x1-5:x1+5]
        noise = abs(np.mean(noise_window))
        StN = cent_counts/noise
        if StN > highest:
            highest = StN
            zero_shift_image = center[i]
            numb = i
            header = fits.getheader(datadir+center[i])

    ## find all shifts for other images:
    imshifts = {} # dictionary to hold the x and y shift pairs for each image
    for image in object_list:
    ## register_translation is a function that calculates shifts by comparing 2-D arrays
        result, error, diffphase = phase_cross_correlation(
            sky_flat_darkcor_data_out[zero_shift_image],
            sky_flat_darkcor_data_out[image])
        imshifts[image] = result



    ## new dictionary for shifted image data:
    shifted_target_data = {}
    for im in object_list:
        ## interp.shift is the function doing the heavy lifting here,
        ## it's reinterpolating each array into the new, shifted one
        shifted_target_data[im] = interp.shift(
            sky_flat_darkcor_data_out[im],
            imshifts[im])


    ## array of aligned arrays:
    targetcube_shift  = np.stack(shifted_target_data.values(),axis=0)

    ## average combined final image:
    target_stacked = np.median(targetcube_shift, axis=0)
    return target_stacked
