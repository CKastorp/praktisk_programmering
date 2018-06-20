.cCFLAGS+=-Wall $$(gsl-config --cflags)
LDLIBS+=$$(gsl-config --libs)

.PHONEY: all

all: plot.svg
plot.svg: plot.gpi output.txt
	gnuplot $<

output.txt:main
	./$< > $@

main: main.o rungekutta.o

.PHONEY: clean
clean:
	-rm -f *.o main output.txt