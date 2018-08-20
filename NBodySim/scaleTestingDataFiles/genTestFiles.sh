#!/bin/bash

dSize=128
while [ $dSize -ge 16 ]; do
  python ./trimTimestep.py $dSize ./untrimmedData $(echo $dSize\startData)
  echo "$dSize complete"
  dSize=$(($dSize-$dSize/4))
done
