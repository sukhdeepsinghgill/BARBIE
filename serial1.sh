#!/bin/bash
#SBATCH -J zbs101                
#SBATCH -p shared                  
#SBATCH -o %j101-out
#SBATCH -e %j101-err
#SBATCH -n 1                      
#SBATCH -t 72:00:00               

LD_LIBRARY_PATH=/home/20ph91f04/local_install/lib/
export LD_LIBRARY_PATH 

exe=./barbie 

$exe 384 101 127 125 1      
