#!/bin/bash

python run.py $1

echo "finished running bcdphot on $1:
`tail -n 50 nohup.out`" | \
mail jlivings@jpl.nasa.gov -s "nohup job finished"
