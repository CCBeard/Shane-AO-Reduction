'''
Flat Fielding
'''

from astropy.io import fits
import numpy as np



def make_master_flats(flat_list, filter, darkcor_data_out, datadir):
    flat_filt_list = []
    for flat in flat_list:
        try:
            flatheader = fits.getheader(datadir+flat)
            flatfilt = flatheader['FILT1NAM']

            if filter == flatfilt:
                flat_filt_list.append(flat) #this makes a list of all flats with this particular filter
            else:
                continue
        except OSError:
            continue

    #making our master flat for this particular filter
    flatcube = np.stack([darkcor_data_out[flat_frame] for flat_frame in flat_filt_list],axis=0)
    master_flat = np.average(flatcube, axis=0)
    normalized_master_flat = master_flat/np.mean(master_flat)

    return normalized_master_flat


def flat_field(object_list, master_flat, darkcor_data_out):
    flat_darkcor_data_out = {}

    for im in object_list:
        #divide each pixel by the flat field
        flat_darkcor_data_out[im] = darkcor_data_out[im]/master_flat

    return flat_darkcor_data_out
