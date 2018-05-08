CFLAGS += $$(gsl-config --cflags)
LDLIBS +=$$(gsl-config --libs)

.PHONEY: all

all: interpolation.svg integration.svg

interpolation.svg: interpolation.gpi output.txt
	gnuplot $<
integration.svg: integration.gpi output.txt
	gnuplot $<
output.txt: main
	./$< > $@

main: main.o linear_interpolation.o findindex.o lininterp_integrator.o

.PHONEY: clean
clean:
	rm -f *.o main output.txt
