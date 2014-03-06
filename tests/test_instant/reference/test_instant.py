#!/usr/bin/env python
#encoding:utf-8

import numpy as np
import netCDF4


def compare_variable(nc1, nc2, var):
    
    var1 = nc1.variables[var]
    var2 = nc2.variables[var]
    
    if var1.shape != var2.shape:
        print 'variable ', var, ' not of the same shape'
        print var1.shape, var2.shape
        return False
    if not np.array_equal(var1[:], var2[:]):
        print 'contents of ', var, ' not equal in both files'
        print 'Average difference : ', np.mean(var1[:]-var2[:])
        print 'Maximum difference : ', np.max(np.abs(var1[:]-var2[:]))
        close = np.isclose(var1[:], var2[:])
        print 'Percentage of similar points : ', 100.*np.sum(close)/np.size(var1)
        return False
        
    print var, ': ok.'
    return True


def test_nonzero(nc2, var):
    
    var2 = nc2.variables[var]

    if np.all(var2[:]==0):
        print 'All points in ', var, ' are zero'
        return False
    print var, ': content appears ok.'

    return True


def compare_instant(file2, file1='./instant_SR_CR_DR_2007-01-01T00-22-49ZN_night_CFMIP2_2.65.nc'):
    
    variables = ['longitude', 'latitude', 'alt_mid', 'alt_bound', 'time', 'SE', 
                'instant_SR', 'instant_CR', 'instant_DR',
                'ATB', 'ATB_mol', 'ATB_per', 'ATB_par', 'TEMP']
        

    print 'Comparing :'
    print file1
    print file2
        
    nc1 = netCDF4.Dataset(file1)
    nc2 = netCDF4.Dataset(file2)
    
    for var in variables:
        compare_variable(nc1, nc2, var)
        
    nc1.close()
    nc2.close()


def test_instant(file2):

    variables = ['longitude', 'latitude', 'alt_mid', 'alt_bound', 'time', 'SE', 
                'instant_SR', 'instant_CR', 'instant_DR',
                'ATB', 'ATB_mol', 'ATB_per', 'ATB_par', 'TEMP']

    print 'Testing content validity of ', file2
    nc2 = netCDF4.Dataset(file2)

    for var in variables:
        test_nonzero(nc2, var)
        
    nc2.close()

def main():

    # file to compare
    file2 = './instant_SR_CR_DR_2007-01-01T00-22-49ZN_night_CFMIP2_2.70.nc'

    compare_instant(file2)
    test_instant(file2)

    print 'tests over.'

if __name__ == '__main__':
    import plac
    plac.call(main)