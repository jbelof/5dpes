
CC=gcc
CFLAGS=-O3
DEFINES=-DDEBUG

.c.o:
	$(CC) -c $(CFLAGS) $(DEFINES) -I. $*.c

all:	main.o cleanup.o input.o pairs.o pbc.o surface.o energy.o
	$(CC) $(CFLAGS) $(DEFINES) *.o -o 5dpes

clean:
	rm -f *.o 5dpes


