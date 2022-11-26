#!/bin/bash
id=2
mkdir -p figs
mkdir -p logs

for m in Batch2 ReMeasure 
do
  for n in 100 
  do
    for n1 in 5 10 15 20 25 30 35 40 45 50
    do 
      for r1 in 1
      do
        for r2 in 0.6
        do
          for a0 in 0 0.25 0.5 
          do
            for a1 in 0.5 
            do
              for e in Gamma t
              do 
                nohup Rscript ./SimNonGaussianS14.R $m $n $n1 $r1 $r2 $a0 $a1 $id $e > ./logs/S1-$m-$n-$n1-$r1-$r2-$a0-$a1-$id-$e.log < /dev/null &
                wait
              done
            done
          done
        done
      done
    done
  done
done