#!/bin/bash

if [ "$#" -ne 2 ]; then
      echo "Creates one list of CALIOP L1 files per day for a given month"
      echo "Usage: $0 YYYY MM"
      exit 1
fi


for i in $(seq -f "%02g" 1 31)
do 
    ./day_list.sh $1 $2 $i
done
