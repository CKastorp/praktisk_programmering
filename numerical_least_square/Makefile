CFLAGS += $$(gsl-config --cflags)
LDLIBS +=$$(gsl-config --libs)

.PHONEY: all

all: fit.svg

fit.svg: fit.gpi output.txt
	gnuplot $<
	
output.txt: main
	./$< > $@

main: main.o QRfactor.o support.o QRsolve.o set_A.o QRleastsquare.o

.PHONEY: clean
clean:
	rm -f *.o main output.txt
