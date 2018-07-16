#!/bin/bash

drawWidth=$(tput cols)
drawHeight=$(tput lines)

let halfWidth=drawWidth/2

if [ $halfWidth -lt $drawHeight ]; then # the window is too narrow
	echo 1
	let drawHeight=$halfWidth
elif [ $halfWidth -gt $drawHeight ]; then # the window is too thicc
	echo "2 $halfWidth $drawHeight"
	let drawWidth=drawHeight*2
fi

echo "Will draw to conole at ${drawWidth}x$drawHeight character resolution"

args="-lm -std=gnu11 -DCONSOLE_W=$drawWidth -DCONSOLE_H=$drawHeight"

if [ $debug = "true" ]; then
	echo "Compiling debug version"
	args="$args -g -O0 -DDEBUG"
else
	printf "Compiling runtime version.\nTo compile for debug, run: 'export debug=true' then then recompile\n"
	args+=" -O3"
fi

if [ `uname` = "Darwin" ]; then
	args+=" -I /opt/X11/include/cairo -L /usr/lib -l cairo"
else
	args+=" -I /usr/include/cairo -L /usr/lib -l cairo"
fi

if [ -z "${debug+x}" ]; then debug="false"; fi

command="gcc ./main.c ./guiOutput.c -o ./simulator $args"
echo "Full compilation instruction is: $command"
eval "$command"

printf "Compilation complete\n"
