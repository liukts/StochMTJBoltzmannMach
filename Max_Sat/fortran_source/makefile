FC           = gfortran
FCFLAGS      = -O3 -Wall -Wextra
F2PYC        = python3 -m numpy.f2py
F2PYFLAGS    = --fcompiler=gnu95 -L. -I.
F2PY_TARGET  = single_sample

#we'll compile the f2py target seperately (using f2py not gfortran)
SOURCES = $(filter-out ./${F2PY_TARGET}.f90, $(wildcard ./*.f90))
OBJS    = $(SOURCES:.f90=.o)

#all, one rule to rule them all
all: ${OBJS}
	${F2PYC} -c ${F2PYFLAGS} ${OBJS} -m ${F2PY_TARGET} ${F2PY_TARGET}.f90
	#in case you want to run the pseudo rngs directly from python as well
	#${F2PYC} -c --fcompiler=gnu95 -m pseudo_rngs pseudo_RNGs.f90

%.o: %.f90
	# $@ evalues to .o
	# $< evaluates to successive .f90
	${FC} ${FCLAGS} -c $< -o $@

# Clean target
clean:
	rm -f ./*.o
	rm -f ./*.mod
	rm -f ./*.so

# Phony targets
.PHONY: all clean
