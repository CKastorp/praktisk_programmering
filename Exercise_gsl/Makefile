CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

.PHONEY: all
all: out.txt airy.svg
	cat $<

main: main.o airy.o

airy.svg: airy.data airy.gpi
	gnuplot airy.gpi

out.txt: main
	./main > $@

clean:
	-rm -f airy.o main.o main airy.data
