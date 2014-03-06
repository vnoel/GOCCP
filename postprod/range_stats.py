#!/usr/bin/env python
#encoding:utf-8

"""
VNoel CNRS 06/03/2014
"""

import glob, os
from subprocess import call

goccp_basedir = '/bdd/CFMIP/CFMIP_OBS_LOCAL/GOCCP'
var_folders = { 'SR_histo330m_':'SR_histo', 
                'SR_histo_Phase330m_':'SR_histo',
                '3D_CloudFraction330m_':'3D_CloudFraction',
                '3D_CloudFraction_Phase330m_':'3D_CloudFraction',
                '3D_CloudFraction_Temp330m_':'3D_CloudFraction',
                'MapLowMidHigh330m_':'MapLowMidHigh',
                'MapLowMidHigh_Phase330m_':'MapLowMidHigh'}
grid_folders = {'CFMIP2':'grid_2x2xL40'}


def move_var(varbase, dir, year, month, dayflag, grid):
    
    print '### Processing ', varbase
    
    files = glob.glob(dir + varbase + '*70.nc')
    if len(files) != 1:
        print '***' + varbase + ' : Cannot find unique file'
        return
    file = files[0]
    print file
    
    target = '/'.join([goccp_basedir, var_folders[varbase], grid_folders[grid], year, dayflag])
    
    rest = file.lstrip(varbase)
    ids = rest.split('_')
    newfilename = varbase + year + '_'.join([month, dayflag, grid, ids[-2], ids[-1]])
    print 'will be moved to ', target + '/' + newfilename

    if not os.path.isdir(target):
        print 'Creating ' + target
        os.makedirs(target)
    
    command = 'cp ' + file + ' ' + target + '/' + newfilename
    print command
    
    call(command, shell=True)
    
    

def main(year='2013', month='01', dayflag='night', grid='CFMIP2'):
    
    period_id = year + month + '_' + dayflag
    print period_id
    
    dirs = glob.glob('../' + period_id + '.*')
    if len(dirs) != 1:
        print '*** Cannot find unique folder'

    dir = dirs[0] + '/out/'
    for var in var_folders:
        move_var(var, dir, year, month, dayflag, grid)


if __name__=='__main__':
    import plac
    plac.call(main)
