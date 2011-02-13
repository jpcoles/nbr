CC=gcc44
CFLAGS=-Wall -O2 -fopenmp
LDFLAGS=-lpng -lm -lgomp

all: nbr

nbr: nbr.o io.o ic.o analysis.o
	$(CC) $(LDFLAGS) nbr.o io.o ic.o analysis.o -o nbr

nbr.o: nbr.c nbr.h
	$(CC) $(CFLAGS) -c nbr.c -o nbr.o

io.o: io.c io.h
	$(CC) $(CFLAGS) -c io.c -o io.o

ic.o: ic.c ic.h
	$(CC) $(CFLAGS) -c ic.c -o ic.o

analysis.o: analysis.c analysis.h
	$(CC) $(CFLAGS) -c analysis.c -o analysis.o


clean: 
	rm -f *.o nbr
