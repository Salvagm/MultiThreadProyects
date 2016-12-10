COMP	:= mpicc
OPT     = -O3 -w -lm

all: mandelbrot_mpi

mandelbrot_mpi: mandelbrot_mpi.c
	$(COMP) $(OPT) -o $@ $@.c
