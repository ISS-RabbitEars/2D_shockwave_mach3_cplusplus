execute:

./run_me.exe

then cd to points dir
then execute the following:

./rename
./render
ffmpeg -framerate 30 -i frame_%d.png sw_30fr.gif
