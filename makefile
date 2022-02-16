.ONESHELL:

SHELL=/bin/bash

CC=g++

CFLAGS=-std=c++2b -fsplit-stack -fno-rtti -fno-exceptions -fsanitize-undefined-trap-on-error -Wall -Wextra -Wdouble-promotion -Wstrict-aliasing=1 -Warith-conversion -Wpointer-arith -Wpedantic -Wno-variadic-macros -Wconversion -Wshadow

INCLUDES=-isystem external/include

DEPS= isingsim.cpp makefile

OPTOPTIONS=-O3
DBGOPTIONS=-g -fsanitize=address -static-libasan

all: isingsim makefile

isingsim: makefile isingsim.cpp external/fmt external/memo
	$(CC) $(CFLAGS) $(INCLUDES) $(OPTOPTIONS) isingsim.cpp -o isingsim

run: isingsim
	rm *.ppm
	./isingsim
	ffmpeg -loglevel warning -pattern_type glob -i "*.ppm" -c:v libx265 -crf 8 -movflags faststart -y output.mov

external/fmt: external/get-fmt.sh
	pushd external > /dev/null
	./get-fmt.sh
	popd > /dev/null

external/memo: external/get-memo.sh
	pushd external > /dev/null
	./get-memo.sh
	popd > /dev/null

clean:
	rm -rf isingsim

debug: $(DEPS) 
	$(CC) $(CFLAGS) $(INCLUDES) $(DBGOPTIONS) isingsim.cpp -o isingsim

