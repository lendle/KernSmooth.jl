WORD_SIZE := $(shell julia -e "print(WORD_SIZE)")

ifeq ($(WORD_SIZE), 64)
	INTFLAG = -fdefault-integer-8
else
	INTFLAG =
endif

OBJS = dqrsl.o dqrdc.o dgefa.o dgesl.o sdiag.o sstdiag.o blkest.o cp.o linbin.o locpoly.o rlbin.o dgedi.o

%.o: %.f Makefile
	gfortran -m$(WORD_SIZE) $(INTFLAG) -fpic -O3 -pipe -g -c $< -o $@

libkernsmooth.so: locpoly.o $(OBJS)
	gcc -std=gnu99 -shared -o libkernsmooth.so $(OBJS) -lblas -lgfortran

clean:
	rm -rf *o
