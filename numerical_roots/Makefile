CFLAGS+=-Wall $$(gsl-config --cflags)
LDLIBS+=$$(gsl-config --libs)

.PHONEY: default

default: output.txt
	cat $<

output.txt:main
	./$< > $@
main: main.o newton.o  QRfactor.o QRsolve.o support.o 

.PHONEY: clean

clean:
	-rm -f *.o main output.txt