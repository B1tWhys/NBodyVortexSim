#!/bin/sh

export debug=false
./compile.sh &&

rm ./data/outputImages/*.png
rm ./data/output.mp4

time ./data/simulator

ffmpeg -f image2 -r 24 -pattern_type glob -i './data/outputImages/*.png' ./data/output.mp4

if [ ${autoplay:-false} = true ]; then
	if [ $(uname) = "Darwin" ]; then
		open ./data/output.mp4
	else
		killall vlc
		vlc -q ./data/output.mp4
	fi
fi