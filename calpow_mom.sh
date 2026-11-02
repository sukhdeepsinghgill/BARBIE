#!/bin/bash

LD_LIBRARY_PATH=/home/20ph91f04/local_install/lib
export LD_LIBRARY_PATH

make calpow_mom


for ((i=101; i<=150; i++))
do
  ./calpow_mom 384 $i
done
