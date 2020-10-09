#!/bin/bash

for n in 200 400 600
do
  for m in 200 400 600
  do
    for d in 2 3
    do
      for dist in poisson binomial
      do
        echo "nohup R CMD BATCH simulation.ants.R $n $m $d $dist &"
        nohup Rscript simulation.ants.R $n $m $d $dist &
      done
    done
  done
  # End of inner loop.
done          
