# -*- MakeFile -*-

# Feb 7, 2022

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
CC = g++
CFLAGS  = -g -Wall -fopenmp -O3 --std=c++11
# the build target(s) executable:
# TARGET = difflin
# If you execute make without a flag, it does the actions under the "all" target by default
# all: difflin32.o
# 	$(CC) $(CFLAGS) -o $(TARGET) difflin32.c simeck32.h
difflin32: difflin32.c difflin32.h
	$(CC) $(CFLAGS) -o difflin32 difflin32.c simeck32.h
difflin48: difflin48.c difflin48.h
	$(CC) $(CFLAGS) -o difflin48 difflin48.c simeck48.h
difflin64: difflin64.c difflin64.h
	$(CC) $(CFLAGS) -o difflin64 difflin64.c simeck64.h
clean:
	rm -f *.o $(TARGET) difflin32 difflin48 difflin64 result*
