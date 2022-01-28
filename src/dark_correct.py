'''
Dark Correction
'''

from astropy.io import fits
import numpy as np


def generate_master_darks(dark_list, drkhdr_list, datadir,):
    '''
    The purpose of this function is to generate a dictionary of master darks, one for each exposure time.
    Usually, when collecting the data, we get several darks for each exposure time we used, unless it was
    a particularly long exposure time. We will then use this master dark to correct all images with the
    same exposure time.
    '''
    darktime_list = [] #list of all the exposure times. There should be at least one for every exposure image we take
    darkname_dict = {}#this will be a dictionary with keys equal to all different exposuretimes. Each entry
    #will list all the dark frames that have that exposure time
    masterdark_dict = {} #this will have one entry for each exposure time: the master dark of that time

    #this populates the darktime list
    for i in range(len(dark_list)):
        try:
            if len(darktime_list) == 0:
                darktime_list.append(drkhdr_list[i])
            elif int(drkhdr_list[i]) not in darktime_list:
                darktime_list.append(drkhdr_list[i])

            else:
                continue
        except OSError:
            continue


    #Now I add the appropriate dark data to each list
    for dark in dark_list:
        try:
            darkheader = fits.getheader(datadir+dark)
            #print(dark)
            darktime = int(float(darkheader['ITIME'])/1000.)
            #print(darktime)
            darkname_dict[darktime].append(dark)
        except KeyError:
            darkname_dict[darktime] = [dark]
        except OSError:
            continue

    #now we go through all the exposure times and make a master dark for each
    for exptime in darktime_list:
        tempdict = {}
        for dark in darkname_dict[exptime]:
            tempdict[dark] = fits.getdata(datadir+dark)
        darkcube = np.stack([tempdict[dark] for dark in darkname_dict[exptime]],axis=0)
        masterdark_dict[exptime] = np.average(darkcube, axis=0)

    return masterdark_dict


def dark_correct(image_name, raw_data, masterdark_dict, datadir):

    header = fits.getheader(datadir+image_name)
    exptime = int(float(header['ITIME'])/1000.)
    obj = header['OBJECT']

    print("Now correcting "+image_name+' using exposure time '+str(exptime)+' s')
    out = raw_data[image_name] - masterdark_dict[exptime]


    #If this is a flat, normalize by exposure time
    if 'flat' in obj.lower():
        out = out/exptime

    return out
