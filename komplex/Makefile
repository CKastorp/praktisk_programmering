CFLAGS = -std=gnu99 -Wall -O
LDLIBS= -lm

.PHONEY: default

default: out.txt	;cat $<

out.txt: main
	./$< > $@
main: main.o komplex.o
komplex.o main.o: komplex.h
clean:
	$(RM) main.o main komplex.o
