#!/bin/bash

python run.py $1

echo "finished running bcdphot on $1:
`cat nohup.out | grep max`" | \
mail jlivings@jpl.nasa.gov -s "nohup job finished"
