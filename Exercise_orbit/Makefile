CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONEY: all
all: logistic_plot.svg orbit.svg

orbit.svg: orbit.gpi output.txt
	gnuplot $<
logistic_plot.svg: logistic_plot.gpi output.txt
	gnuplot $<
output.txt: main
	./$< > $@
main: main.o logistic.o relativistic.o

.PHONEY: clean
clean:
	$(RM) logistic.o main.o main output.txt logistic_plot.svg relativistic.o orbit.svg
