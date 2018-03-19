CFLAGS += -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

.PHONEY: testplot.svg
testplot.svg: testplot.gpi cosines.txt sines.txt
	gnuplot testplot.gpi

cosines.txt: main-stdio
	echo 1 2 0.5 4 1.2 2.3 |./$< > $@
main-stdio: main-stdio.o

main-stdio.o: main-stdio.c

sines.txt: main-cmdline
	./$< 0 1.2 0.2 2.3 0.7 1.4 3 > $@

main-cmdline: main-cmdline.o

main-cmdline.o: main-cmdline.c
clean:
	$(RM) main-stdio.o main-cmdline.o cosines.txt sines.txt testplot.svg
