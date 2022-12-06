#!/bin/bash
for d in 100 200 400
    do
for n in 5 10 25 
    do
for T in 2 4 6 
    do         
for ord in 0 1 2 4
    do
    sbatch --export arg1=$d,arg2=$n,arg3=$T,arg4=$ord thm.cmd
done
done
done 
done