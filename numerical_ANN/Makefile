CFLAGS+=-Wall $$(gsl-config --cflags)
LDLIBS+=$$(gsl-config --libs)
CC=gcc

.PHONEY: all

all: ANN.svg
ANN.svg: network.gpi output.txt
	gnuplot $<

output.txt:main
	./$< > $@

main: main.o ANN.o ANN.h min.o

.PHONEY: clean
clean:
	-rm -f *.o main output.txt