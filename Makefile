# Compiler and flags
FC=ifort
FCFLAGS=-c -qopenmp #-qopt-report=5
FLFLAGS=-qopenmp

# source files and objects
SRCS = $(patsubst %.f, %.o, $(wildcard *.f))

# Program name
PROGRAM = ldgc

all: $(PROGRAM)

run:
		@echo "******************* Running ***********************"
		time ./$(PROGRAM)

runB: $(PROGRAM)
		@echo "******************* Running ***********************"
		time ./$(PROGRAM)

# Linking
$(PROGRAM): $(SRCS)
		@echo "******************* Linking ***********************"
		$(FC) $(FLFLAGS) -o $@ $<
		@echo "Linking: OK"

# Source compilation
%.o: %.f
		@echo "******************* Compilation *******************"
	    $(FC) $(FCFLAGS) -o $@ $<
		@echo "Compilation: OK"

clean:
		@echo "******************* Cleaning **********************"
		rm -f *.optrpt *.o ldgc *~ fort* *.eco *.con *.optrpt
		@echo "Removing: OK"
