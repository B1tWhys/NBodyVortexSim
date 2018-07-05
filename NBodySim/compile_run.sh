#!/bin/sh

export debug=false
./compile.sh &&

rm ./outputImages/*.png
rm ./output.mp4

time ./simulator

ffmpeg -f image2 -r 24 -pattern_type glob -i './outputImages/*.png' output.mp4

# if [ $(uname) = "Darwin" ]; then
# 	open ./output.mp4
# else
# 	killall vlc
# 	vlc -q ./output.mp4
# fi