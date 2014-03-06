#!/usr/bin/env python
#encoding:utf-8

"""
VNoel CNRS 06/03/2014
"""

import glob, os
from subprocess import call

goccp_basedir = '/bdd/CFMIP/CFMIP_OBS_LOCAL/GOCCP'
grid_folders = {'CFMIP2':'grid_L40'}


def move_instant(dir, year, month, grid):
    
    mask = dir + '*.nc'
    files = glob.glob(mask)
    if len(files) == 0:
        print '*** Cannot find files'
        return

    folder = year + month
    target = '/'.join([goccp_basedir, 'instant_SR_CR_DR', grid_folders[grid], year, folder])

    print '%d instant files will be moved in ' % len(files), target
    if not os.path.isdir(target):
        print 'Creating ' + target
        os.makedirs(target)
    
    command = 'cp ' + mask + ' ' + target + '/'
    print command
    
    call(command, shell=True)
    
    

def main(year='2013', month='01', dayflag='night', grid='CFMIP2'):
    
    period_id = year + month + '_' + dayflag
    print period_id
    
    dirs = glob.glob('../' + period_id + '.*')
    if len(dirs) != 1:
        print '*** Cannot find unique folder'

    dir = dirs[0] + '/out/instant/'
    move_instant(dir, year, month, grid)


if __name__=='__main__':
    import plac
    plac.call(main)
