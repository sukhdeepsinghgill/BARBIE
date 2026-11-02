#!/bin/bash

for i in {101..150}; do
  sbatch serial$i.sh
done
