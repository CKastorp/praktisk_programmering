CFLAGS =-std=gnu99 -Wall -o
LDLIBS = -lm

.PHONEY: all
all: out.txt
	cat out.txt
out.txt: main
	 ./main > out.txt
#math: math.o
#	cc math.o -o math
main: main.o
	$(CFLAGS) cc main.o -o main
main.o: main.c
	$(CFLAGS) cc -c main.c $(LDLIBS)-o main.o
clear:
	$(RM) main.o out.txt
