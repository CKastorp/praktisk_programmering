CFLAGS += -Wall -std=gnu99
CFLAGS += $$(gsl-config --cflags)
LDLIBS += $$(gsl-config --libs)

all: output ; cat $<
output: main ; ./$< > $@ $(CFLAGS) $(LDLIBS)
main: main.o
main.o: main.c
clean: ; $(RM) main main.o output
