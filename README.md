GOCCP
=====

source code for GOCCP production - LMD branch

# 1. How to compile

./compile.sh calmdz2.70 on 0

* calmdz2.70 = name of executable.
* on/off = produce instant files
* 0/1 = non-overlap mode if 1

This will compile the source and link the executable in the run/ directory.

# 2. How to run

## 2.1 On a given list of files

* enter the run directory
* create in IN/ an input file listing the CALIOP L1 files to process
* ./calmdz2.70
* enter requested arguments
* output will be in ./out

## 2.2 On a given list of files

* create in run/IN an input file listing the CALIOP L1 files to process
* qsub run_list.pbs -v liste=input_list_file
* output will be in input_list_file.$JOBID/out/

## 2.3 For day segments over a month

* create daily lists using preprod/day_lists.sh YYYY MM (they go in run/IN)
* ./run_month.sh YYYY MM DAYFLAG (DAYFLAG is day or night)
* output will be in run.YYYYMM/out/
