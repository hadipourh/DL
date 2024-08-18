# -*- MakeFile -*-
# July 25, 2023

# compiler
# define the C compiler to use
# for C++ define  CC = g++
CC = g++
CFLAGS  = -g -Wall -fopenmp -O3

all: difflin differential linear

difflin: difflin.c
	$(CC) $(CFLAGS) -o difflin difflin.c

differential: differential.c
	$(CC) $(CFLAGS) -o differential differential.c

linear: linear.c
	$(CC) $(CFLAGS) -o linear linear.c

clean:
	rm -f difflin differential linear

