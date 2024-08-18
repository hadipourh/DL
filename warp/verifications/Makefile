# -*- MakeFile -*-
# Feb 9, 2022

# compiler
# define the C compiler to use
# for C++ define  CC = g++
CC = g++
CFLAGS  = -g -Wall -fopenmp -O3 --std=c++11


all: difflin diff

difflin: difflin.c warp.c warp.h difflin.h
	$(CC) $(CFLAGS) -o difflin difflin.c warp.c

diff: diff.c warp.c warp.h
	$(CC) $(CFLAGS) -o diff diff.c warp.c
clean:
	rm -f difflin diff result*
