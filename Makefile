CC=gcc -fopenmp
CFLAGS= -g -Ofast -march=native -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=leak\
         `pkg-config --cflags libconfig`
#maybe delete -fsanitize=address
#CFLAGS = -O2 `pkg-config --cflags libconfig` `gsl-config --cflags`
DEBUG = -Wall -Wextra -pedantic -ggdb -D_DEBUG_
INCLUDE = include

COMPILE = $(CC) $(CFLAGS) $(DEBUG) -I$(INCLUDE)
LIBS= -lfftw3 -lm `pkg-config --libs libconfig` #lfftw3 is needed for FFT

src = $(wildcard src/*.c)
obj = $(src:.c=.o)
dep = $(obj:.o=.d)

all: vortices

vortices: $(obj)
	$(COMPILE) $^ -o $@ $(LIBS)

-include $(dep)

%.d: %.c
	$(CPP) $(CFLAGS) $< -I$(INCLUDE) -MM -MT $(@:.d=.o) >$@

%.o: %.c
	$(COMPILE) -c $< -o $@

.PHONY: clean
clean:
	rm vortices $(obj) $(dep)
