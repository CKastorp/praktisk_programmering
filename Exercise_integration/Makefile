CFLAGS+= -Wall $$(gsl-config --cflags)
LDLIBS =$$(gsl-config --libs)

.PHONEY: all
all:energies.svg
energies.svg: energies.gpi output.txt
	gnuplot $<
output.txt: main
	./$< > $@
main: main.o integration.o dirac.o

.PHONEY: clean
clean:
	$(RM) energies.svg main.o integration.o dirac.o main output.txt
