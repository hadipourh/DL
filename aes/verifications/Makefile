all:
	gcc aesni.h difflin.c -g -O1 -lm -finline-functions -fomit-frame-pointer -funroll-loops -Wall -msse2 -msse -march=native -maes -o  difflin
clean:
	rm -f *.o difflin
