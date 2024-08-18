# -*- MakeFile -*-

# @H. Hadipour
# Sep 19, 2023

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

# Compiler and flags
CC = gcc
CFLAGS = -g -Wall

# Target executable
TARGET = serpentopt

# Default target
all: $(TARGET)

# Linking rule
$(TARGET): serpentopt.o
	$(CC) $(CFLAGS) -o $(TARGET) serpentopt.o

# Compilation rule
serpentopt.o: serpentopt.c
	$(CC) $(CFLAGS) -c serpentopt.c

# Clean rule
clean:
	rm -f *.o $(TARGET)