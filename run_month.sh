#!/bin/bash

if [ "$#" -ne 3 ]; then
      echo "Runs GOCCP on the cluster on a full month, one CPU per day."
      echo "List files for each day must exist (cf preprod/day_lists.sh)"
      echo "Usage: $0 YYYY MM DAYFLAG"
      exit 1
fi

YEAR=$1
MONTH=$2
dayflag=$3

workdir=run.$YEAR$MONTH
mkdir $workdir
cd $workdir

ln -s ../run/* .
ln -s ../preprod/run_month.pbs .
mkdir -p out/instant

qsub -v year=$YEAR,month=$MONTH,dayflag=$dayflag run_month.pbs