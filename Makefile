# Makefile for C++ files

compile:
	g++ -std=c++17 -O3 -o code startup_code.cpp

clean:
	rm -f code

run:
	./code

all:
	make clean
	make compile
	make run

gdb:
	g++ -std=c++17 -g -Wall -o code startup_code.cpp
	gdb code