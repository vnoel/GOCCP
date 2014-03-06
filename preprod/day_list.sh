#!/bin/bash

if [ "$#" -ne 3 ]; then
      echo "Creates a list of CALIOP L1 files for a given date"
      echo "Usage: $0 YYYY MM DD"
      exit 1
fi

find /bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.02/$1/$1_$2_$3* | grep ZN.hdf > IN/$1$2$3_night
