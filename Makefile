CC=g++
#CFLAGS=-Wall -c
CFLAGS=-Wall -c -std=c++11
CFLAGS1=-Wall -c -std=c++11

all: check


check: check.o Matrix.o isomorphi.o
	$(CC) check.o Matrix.o isomorphi.o -o check

check.o: check.cpp
	$(CC) $(CFLAGS) check.cpp

isomorphi.o: isomorphi.cpp
	$(CC) $(CFLAGS1) isomorphi.cpp

#semigroups: semigroups.o Matrix.o
#	$(CC) Matrix.o semigroups.o -o semigroups
#
#semigroups.o: semigroups.cpp
#	$(CC) $(CFLAGS) semigroups.cpp
#
#graphs: graphs2.o Matrix.o
#	$(CC) Matrix.o graphs2.o -o graphs
#
#graphs2.o: graphs2.cpp
#	$(CC) $(CFLAGS) graphs2.cpp

Matrix.o: Matrix.cpp
	$(CC) $(CFLAGS) Matrix.cpp

clean:
	rm -rf *.o check
