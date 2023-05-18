'''
Sky Subtraction
'''

from astropy.io import fits
import numpy as np
from utils import guess_gaussian_parameters, weighted_centroid


def make_master_sky(object_list, flat_darkcor_sigmacut_data, datadir):
    ## create an empty dictionary to populate with the completely corrected science frames:
    cr_flat_darkcor_data_out = {}

    # let's add some bits to account for different exposure times on images
    exposuretimes = []  # this will list all the different exposure times used on this image
    for obj in object_list:
        header = fits.getheader(datadir + obj)
        if len(exposuretimes) == 0:
            exposuretimes.append(float(header['ITIME']) / 1000)
        elif float(header['ITIME']) / 1000 not in exposuretimes:
            exposuretimes.append(float(header['ITIME']) / 1000)
        else:
            continue
    exp_dict = {}
    # make a list of each data image in each exposure time
    # skycube = []
    for time in exposuretimes:
        newlist = []
        for i in range(len(object_list)):
            header = fits.getheader(datadir + object_list[i])
            if time == float(header['ITIME']) / 1000:
                newlist.append(object_list[i])
            else:
                continue
        exp_dict[time] = newlist

    master_sky_dict = {}  # dictionary of master skys for each exposure time
    for time in exposuretimes:
        print('Sky Subtracting images with exp time = ' + str(time) + ' s')
        # this next section is for use with a dither script, and places each object image in one of the 5 dither positions
        btmrt = []
        btmlft = []
        toprt = []
        toplft = []
        center = []
        ### we assume the first image is taken at center position, we use that to constrain the positions.
        first_im = min([int(key[1:5]) for key in flat_darkcor_sigmacut_data.keys()])
        x1, y1 = weighted_centroid(flat_darkcor_sigmacut_data[f's{first_im:04d}.fits'])
        # print('center position: ', x1, y1)
        for im in exp_dict[time]:  # for now just checking via centroid position
            # x0,y0,sigma,A = guess_gaussian_parameters(flat_darkcor_sigmacut_data[im])
            x0, y0 = weighted_centroid(flat_darkcor_sigmacut_data[im])
            # print(x0, y0)
            if x0 > x1 + 10:
                if y0 < y1:
                    btmrt.append(im)
                elif y0 > y1:
                    toprt.append(im)
            elif x0 < x1 - 10:
                if y0 < y1:
                    btmlft.append(im)
                elif y0 > y1:
                    toplft.append(im)
            elif y0 > y1 + 10:
                if x0 < x1:
                    toplft.append(im)
                elif x0 > x1:
                    toprt.append(im)
            elif y0 < y1 - 10:
                if x0 < x1:
                    btmlft.append(im)
                elif x0 > x1:
                    btmrt.append(im)
            else:
                center.append(im)

        positions = [len(btmrt), len(btmlft), len(toprt), len(toplft), len(center)]
        nonzero = []
        pos = 5
        for num in positions:
            if num != 0:
                nonzero.append(num)
            else:
                pos = pos - 1
                continue

        min_dither = np.min(nonzero)  # wanna know the smallest number of nonzero positions

        # If we don't have enough dither positions, this doesn't work
        print('Positions = ' + str(pos))
        if pos < 3:
            print('Exposure Time: ' + str(exposuretimes[i]))
            print('Gonna delete: ')
            print(explist[i])
            for obj in exp_dict[time]:
                object_list.remove(obj)
                continue

        skycubelist = np.concatenate(
            [toprt[0:min_dither], toplft[0:min_dither], btmrt[0:min_dither], btmlft[0:min_dither],
             center[0:min_dither]])
        skycube = np.stack([flat_darkcor_sigmacut_data[science_frame] for science_frame \
                            in skycubelist], axis=0)
        master_sky_dict[time] = np.median(skycube, axis=0)

    return master_sky_dict, exp_dict, center


def sky_subtract(object_list, exptimes, flat_darkcor_sigmacut_data, exp_dict, master_sky_dict):
    sky_flat_darkcor_data_out = {}
    for time in exptimes:
        for im in exp_dict[time]:
            skysub_image = flat_darkcor_sigmacut_data[im] - master_sky_dict[time]
            sky_flat_darkcor_data_out[im] = skysub_image

        return sky_flat_darkcor_data_out
