makeSim: main.c main.h constants.h
ifeq ($(debug), true)
		gcc ./main.c -O0 -g -lm -std=gnu11 -o ./simulator
else
		gcc ./main.c -O3 -lm -std=gnu11 -o ./simulator
endif
