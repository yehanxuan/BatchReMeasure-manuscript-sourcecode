#!/bin/bash

id=2

mkdir -p figs
mkdir -p logs

for m in Ignore Batch2 ReMeasure
do
  for nc1 in 50 100 200
  do
    for nc2 in $((nc1/10)) $((nc1/5)) $((nc1/2))
    do 
      for s1 in 1
      do
        for s2 in 2
        do
          for rho in 0.5 0.8
          do
            for a0 in 0 1 2
            do
              for a1 in 0.5
              do
              nohup Rscript ./Main_S1.R $m $nc1 $((nc1/2)) $nc2 $s1 $s2 $rho $a0 $a1 $id > ./logs/S1-$m-$nc1-$nc1-$nc2-$s1-$s2-$rho-$a0-$a1-$id.log < /dev/null &
              wait
              done
            done
          done
        done
      done
    done
  done
done