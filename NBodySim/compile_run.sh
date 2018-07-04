#!/bin/sh

export debug=false
./compile.sh &&

rm ./outputImages/*.png
rm ./outputImages/*.mp4

./simulator

cd ./outputImages

ffmpeg -f image2 -r 24 -pattern_type glob -i '*.png' output.mp4

open ./output.mp4