CFLAGS+=-Wall $$(gsl-config --cflags)
LDLIBS+=$$(gsl-config --libs)

.PHONEY: all

all: output.txt
	cat $<

output.txt:main
	./$< > $@

main: main.o int2D.o

.PHONEY: clean
clean:
	-rm -f *.o main output.txt