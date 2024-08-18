# -*- MakeFile -*-
# July 8, 2022

# compiler
# define the C compiler to use
# for C++ define  CC = g++
CC = g++
CFLAGS  = -g -Wall -fopenmp -O3


all: difflin

lintest: difflin.c
	$(CC) $(CFLAGS) -o difflin difflin.c

clean:
	rm -f difflin
