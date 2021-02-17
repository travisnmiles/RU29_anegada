import numpy as np

def grid_glider(dataset_id,
               varz2d = ['potential_temperature', 'salinity', 'cdom', 'chlorophyll_a', 'beta_700nm'],
                zgrid = np.arange(0,1000,5),
               ):
    '''grid the glider data from RUCOOL Erddap. this needs work'''
    import xarray as xr
    import pandas as pd
    from erddapy import ERDDAP
    
    from scipy.signal import find_peaks
    from scipy import stats
    e = ERDDAP(
    server="http://slocum-data.marine.rutgers.edu/erddap",
    protocol="tabledap",
    response="nc",
    )

    # get the science data:
    e.dataset_id = dataset_id

    # this connects to the data and load into an pandas dataframe
    ds = e.to_pandas()
    # remove the spaces from the column names
    ds.columns = ds.columns.str.split(' ').str[0]

    # get the time to be a datetime object
    ds['time'] = pd.to_datetime(ds['time'])

    # put the times in order
    ds = ds.sort_values(by=['time'])
    
    # fill nans in dpeth for the profile breakup
    interpd = ds.depth.interpolate()
    
    # find the top and bottom of each profile
    apogee, prop = find_peaks(interpd.values,  threshold=None, distance=None, prominence=50)

    perogee, prop = find_peaks(-1*interpd.values,  threshold=None, distance=None, prominence=50)

    # stack the index of the turning points into one vector
    turns = np.sort(np.append(apogee, perogee ))
    
    
    # this is your depth grid, you can set:
    zgrd = zgrid

    # list of variables to grid in 2d:
    # you choose from the columns of the science data
    dataz = varz2d


    # this is a dict to hold our gridded stuff
    # until we make a dataset later
    d2 = {}

    # loop on the variables you want to bin
    for varz in dataz:    
        values = ds[varz] # grab some data

        #this thing below bins the data
        ret = stats.binned_statistic_2d(ds.index.values, ds.depth, values, statistic='mean', bins=[ turns, zgrd ])
        d2[varz] = ret.statistic.T
        
    # things to bin in the x direction
    oneDvars = ['latitude','longitude', 'time', 'u', 'v']

    # NB: u, v only have one value per dive sequence, so only half the number profiles!
    # actually, its weirder than that... not sure there are more than half...

    # dict to hold our 1d bins
    d1 = {}

    # loop on 1d stuff:
    for thing in oneDvars:    
        if thing == 'time':
            bin_means, bin_edges, binnumber = stats.binned_statistic(ds.index.values,
                        ds[thing].astype(int), statistic = 'mean', bins=turns)
            bin_means = pd.to_datetime(bin_means)
        else:

            bin_means, bin_edges, binnumber = stats.binned_statistic(ds.index.values,
                        ds[thing].values, statistic = np.nanmean, bins=turns)
        d1[thing] = bin_means
        
    # need the depth grid centers
    zgrd_ctr = zgrd[:-1] + np.diff(zgrd).mean()/2

    # create the dataset
    ds_gridded = xr.Dataset( coords = {'date': d1['time'].values,'depth': zgrd_ctr ,
                               'lat': ('date', d1['latitude']),
                               'lon': ('date', d1['longitude'])
                              },
                   data_vars = {'u': ('date', d1['u']), 
                               'v': ('date', d1['v'])})

    # add the other data
    for varz in dataz:
        ds_gridded[varz] = ( ('depth', 'date'),d2[varz] )
        
    return ds_gridded




