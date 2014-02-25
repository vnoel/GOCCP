#!/usr/bin/env python

import netCDF4
import os, sys

import matplotlib.pyplot as plt

def maps_show(ncfile, variables):

    print 'Opening ', ncfile

    nc = netCDF4.Dataset(ncfile)
    try:
        lon = nc.variables['longitude'][:]
        lat = nc.variables['latitude'][:]
    except:
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
            
    for var in variables:
        plt.figure()
        data = nc.variables[var][0,...]
        plt.pcolormesh(lon, lat, data)
        plt.colorbar()
        plt.title(var)
        
    plt.show()
    

def main(filepath):
    
    ncfile = os.path.basename(filepath)
    
    if ncfile.startswith('Map'):
        if ncfile.startswith('MapLowMidHigh_Phase'):
            sys.exit('Sorry !')
        elif ncfile.startswith('MapLowMidHigh'):
            variables = ['cllcalipso', 'clmcalipso', 'clhcalipso', 'cltcalipso', 'clccalipso']
        elif ncfile.startswith('MapHigh'):
            variables = ['clhcalipso', 'bascalipso', 'topcalipso']

        maps_show(filepath, variables)


if __name__ == '__main__':
    import plac
    plac.call(main)