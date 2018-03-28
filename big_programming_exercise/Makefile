CFLAGS +=-Wall $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.Phoney: all
all: report.pdf
report.pdf: report.tex plot.pdf
	pdflatex $<
	pdflatex $<
plot.pdf: plot.gpi output.txt
	gnuplot $<
output.txt: main
	./$< > $@
main: main.o

.Phoney: clean

clean:
	$(RM) main.o main output.txt plot.pdf report.pdf report.aux report.log
