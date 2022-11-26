#!/bin/bash

id=2

mkdir -p figs
mkdir -p logs

for m in Batch2 ReMeasure Ignore LS 
do
  for n in 100 
  do
    for n1 in 5 10 15 20 25 30 35 40 45 50
    do 
      for r1 in 2 
      do
        for r2 in 0.3 0.6 0.9  
        do
          for a0 in 0 0.5 0.8  
          do
            for a1 in 0.5 
            do
              nohup Rscript ./Sim3.2.R $m $n $n1 $r1 $r2 $a0 $a1 $id > ./logs/S1-$m-$n-$n1-$r1-$r2-$a0-$a1-$id.log < /dev/null &
              wait
            done
          done
        done
      done
    done
  done
done

