# CFLAGS  = -g -std=c99 -O3

all: outro

outro:
	gcc ray_paralel.c -o ray_p -lm -O3 -fopenmp
	gcc ray_paralel.c -o ray -lm -O3

clean:
	rm ray ray_p image.ppm 
