# -*- MakeFile -*-

# Feb 9, 2022

# The main structure of a make file commands is as follows: 
# target: dependencies
#       action
# Example: 
# output: main.o minor1.o minor2.o
#         g++ -o main main.o minor1.o minor2.o
# Some compiler flags:
# -E                       Preprocess only; do not compile, assemble or link.
# -S                       Compile only; do not assemble or link. (generates an assembly file and stop)
# -c                       Compile and assemble, but do not link. (generates an assembly code and converts it to the machine code by assembler)
# -o <file>                Place the output into <file>. (Links assembled files, and builds an executable file)
# -g    adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
# If you don't use a flag at all, your compiler (gcc here) do the process, compile

# compiler
# define the C compiler to use
# for C++ define  CC = g++
CC = gcc
CFLAGS  = -g -Wall -O3 -lm
# the build target(s) executable:
TARGET = all
# If you execute make without a flag, it does the actions under the "all" target by default
all: difflin

difflin: difflin.c 
	$(CC) -o difflin difflin.c $(CFLAGS)

clean:
	rm -f *.o $(TARGET) difflin