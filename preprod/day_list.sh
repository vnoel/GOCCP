#!/bin/bash

if [ "$#" -ne 4 ]; then
      echo "Creates a list of CALIOP L1 files for a given date"
      echo "Usage: $0 YYYY MM DD ZN/ZD"
      exit 1
fi

if [ $4 == 'ZN' ]; then
    dayflag='night'
else
    dayflag='day'
fi

find /bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.30/$1/$1_$2_$3* | grep $4.hdf > IN/$1$2$3_$dayflag
