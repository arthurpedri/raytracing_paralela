# CFLAGS  = -g -std=c99 -O3

all: outro

outro:
	mpicc ray_mpi.c -o ray -lm -O3

clean:
	$(RM) ray *.o
