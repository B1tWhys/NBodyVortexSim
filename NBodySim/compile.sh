#!/bin/bash

drawWidth=$(tput cols)
drawHeight=$(tput lines)

let doubleWidth=drawWidth*2

if [ $doubleWidth -lt $drawHeight ]; then
	let drawHeight=drawWidth/2
elif [ $doubleWidth -gt $drawHeight ]; then
	let drawWidth=drawHeight*2
fi

echo "Will draw to conole at ${drawWidth}x$drawHeight character resolution"

args="-lm -std=gnu11 -DCONSOLE_W=$drawWidth -DCONSOLE_H=$drawHeight"

if [ `uname` = "Darwin" ]; then
	args+=" -I /opt/X11/include/cairo -L /usr/lib -l cairo"
else
	args+=" -I /usr/include/cairo -L /usr/lib -l cairo"
fi

if [ -z "${debug+x}" ]; then debug="false"; fi

if [ $debug = "true" ]; then
	echo "Compiling debug version"
	args="$args -DDEBUG -g -O0"
else
	printf "Compiling runtime version.\nTo compile for debug, run: 'export debug=true' then then recompile\n"
	# args="$args -O3"
fi

command="gcc ./main.c ./guiOutput.c -o ./simulator $args"
echo "Full compilation instruction is: $command"
eval "$command"

printf "Compilation complete\n\n\n"
