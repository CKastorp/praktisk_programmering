CFLAGS+= -Wall $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.Phoney:all
all: path.svg swave.svg out.txt
	sed -n -e 276,280p -e 1085p  out.txt
swave.svg: swave.gpi out.txt
	gnuplot $<
path.svg: path.gpi out.txt
	gnuplot $<
out.txt: main
	./$< > $@

main: main.o rosenbrock.o hydrogen.o

.Phoney: clean
clean:
	$(RM) main.o rosenbrock.o output.txt hydrogen.o path.svg
