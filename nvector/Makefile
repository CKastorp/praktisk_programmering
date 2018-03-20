.PHONEY: all

all: output.txt
	cat $<
output.txt: main-nvec
	./$< > $@

main-nvec: main-nvec.o Nvector.o
	$(CC) $^ -o $@
main-nvec.o Nvector.o: Nvector.h

.PHONEY: clean
clean:
	$(RM) main-nvec.o Nvector.o main-nvec output.txt
