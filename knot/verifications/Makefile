# -*- MakeFile -*-
# July 25, 2023

# compiler
# define the C compiler to use
# for C++ define  CC = g++
CC = g++
CFLAGS  = -g -Wall -fopenmp -O3


all: difflin

diff-lin: difflin.c
	$(CC) $(CFLAGS) -o difflin difflin.c

diff: diff.c
	$(CC) $(CFLAGS) -o diff diff.c
clean:
	rm -f difflin diff result*.txt
