#!/bin/bash

for inFile in ./*startData; do
  fSize=$(echo $inFile | ack -o "\d*")
  
  echo "domainSize $fSize $(../data/simulator $fSize ./inFile)"`# | tee -a ./testResults`
  printf "$(../data/simulator $fSize)"
done
