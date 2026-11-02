CC = gcc
CFLAGS ?= -O2 -std=gnu11 -Wall -Wextra -fopenmp
LDLIBS = -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm

.PHONY: all clean
all: barbie field calpow_mom

barbie: barbie.c funcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) barbie.c funcs.c $(LDFLAGS) $(LDLIBS) -o $@

field: generate_field.c ngrf_2d_funcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) generate_field.c ngrf_2d_funcs.c $(LDFLAGS) $(LDLIBS) -o $@

calpow_mom: calpow_mom.c funcs.c
	$(CC) $(CPPFLAGS) $(CFLAGS) calpow_mom.c funcs.c $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f barbie field calpow_mom
