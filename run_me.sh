#!/bin/bash

mkdir points
g++ swp.cpp -o test
./test
cp render.sh ./points/
cp points_template.pov ./points/
cd points
./render.sh
ffmpeg -framerate 30 -i frame_%d.png sw_30fr.gif
ffmpeg -framerate 30 -i frame_%d.png sw_30fr.mp4
