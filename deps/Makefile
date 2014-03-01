WORD_SIZE := $(shell julia -e "print(WORD_SIZE)")

ifeq ($(WORD_SIZE), 64)
	INTFLAG = -fdefault-integer-8
else
	INTFLAG =
endif


%.o: %.f Makefile
	gfortran -m$(WORD_SIZE) $(INTFLAG) -fpic -O3 -pipe -g -c $< -o $@

libkernsmooth.so: locpoly.o linbin.o rlbin.o dgefa.o dgesl.o
	gcc -std=gnu99 -shared -o libkernsmooth.so dgefa.o dgesl.o linbin.o locpoly.o rlbin.o -lblas -lgfortran

clean:
	rm -rf *o
