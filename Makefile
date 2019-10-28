all: minibb.o

FLAGS=-Wall -pedantic -std=c++17 -O3
INCLUDES=-I/usr/include/eigen3

minibb.o: curves.cc curves.hh vector.hh
	g++ -c -o $@ $< $(FLAGS) $(INCLUDES)
