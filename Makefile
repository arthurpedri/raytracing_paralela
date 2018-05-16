# CFLAGS  = -g -std=c99 -O3

all: outro

outro:
	gcc ray.c -o ray -lm -O3

clean:
	$(RM) ray *.o
