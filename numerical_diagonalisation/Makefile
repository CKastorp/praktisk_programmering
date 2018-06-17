CFLAGS+= -Wall $$(gsl-config --cflags)
LDLIBS+=$$(gsl-config --libs)
CC=cc

.PHONEY: default

default: output.txt time.svg
	cat $<
output.txt: main
	./$< 6 >$@
main: main.o support.o jacobi_decomp.o

time.svg: CPUtime.gpi time.txt
	gnuplot $<

data time.txt: main Makefile
	cat /dev/null > $@
	for n in  `seq 50 50 1000`; do \time --format "$$n %U" --append --output time.txt ./$< $$n >data; done

.PHONEY: clean

clean:
	-rm -f *.o main output.txt time.txt time.svg