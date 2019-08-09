######################################################
# Authors: Fahimeh Baftizadeh, Taylor Connington
# Date created: 9/1/2017
######################################################

"""
Place to hold analysis functions. 
Any transformation / operation of the data tables would fit.
"""
import numpy as np

#################################################
#
#     Analysis
#
#################################################

def spherical_coords(df):
    rho = np.sqrt(df.x ** 2 + df.y ** 2)
    theta = np.arctan2(rho, df.z)
    phi = np.arctan2(df.y, df.x)

    return theta, phi

def find_thresholds(df):
    thresholds = {}
    g = df.groupby('electrode')

    for el, group in g:

        # find threshold
        sub = group[['num_spikes', 'amp']]

        # get index of threshold edges based on amp
        # using num_spikes values
        spike = sub[sub['num_spikes'] > 0]
        nospike = sub[sub['num_spikes'] == 0]

        # use indexes to obtain amp values, if any
        if spike.size == 0 or nospike.size == 0:
            # print 'cannot find threshold range for el {}'.format(el)
            thresh = None
        else:
            thresh_upper = group.loc[spike.idxmax().amp].amp    # get the most positive/ least neg. amp
            thresh_lower = group.loc[nospike.idxmin().amp].amp  # get the most negative amp
            thresh = [thresh_lower, thresh_upper]

        thresholds[el] = thresh

    return thresholds