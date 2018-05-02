CFLAGS+=-Wall $$(gsl-config --cflags)
LDLIBS+= $$(gsl-config --libs)

.PHONEY: all
all: path.svg fit.svg output.txt
	sed -n '173,176p' output.txt
fit.svg: fit.gpi output.txt
	gnuplot $<

path.svg: path.gpi output.txt
	gnuplot $<
output.txt: main
	./$< > $@
main: main.o findmin.o curvefit.o

.PHONEY: clean
clean:
	$(RM) main *.svg *.o output.txt
