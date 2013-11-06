CFLAGS=-I. -DDEBUG -O3 -g
#CFLAGS=-I. -O3 -march=native -mtune=native
LDLIBS=-lgsl -lgslcblas

test_weierstrass: test_weierstrass.c weierstrass.c weierstrass.h

clean:
	rm -f test_weierstrass *.o
