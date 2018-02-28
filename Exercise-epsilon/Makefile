CFLAGS = -std=gnu99 -Wall -O

.PHONEY: default

default: out.txt
	cat out.txt

out.txt: main
	./main > out.txt
main: maxima.o epsilon.o summation.o main.o equal.o
	$(CC) $(CFLAGS) $^ -o main -lm
main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o
summation.o: summation.c
	$(CC) $(CFLAGS) -c summation.c -o summation.o
equal.o: equal.c
	$(CC) $(CFLAGS) -c equal.c -o equal.o
maxima.o: maxima.c
	$(CC) $(CFLAGS) -c maxima.c -o maxima.o
epsilon.o: epsilon.c
	$(CC) $(CFLAGS) -c epsilon.c -o epsilon.o
clean:
	$(RM) out.txt epsilon.o summation.o main.o
