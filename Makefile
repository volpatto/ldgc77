# Compiler and flags
FC=ifort
#FCFLAGS=-c -openmp -vec-report5 -vec-report-file=opt_report.txt
FCFLAGS=-c -openmp -opt-report-phase=hpo -opt-report-file=opt_report.txt
FLFLAGS=-openmp

#Threads
NUMTHREADS=4
NUMTHREADS?=

# source files and objects
SRCS = $(patsubst %.f, %.o, $(wildcard *.f))

# Program name
PROGRAM = ldgc

all: $(PROGRAM)

run:
		@echo "******************* Running ***********************"
		time ./$(PROGRAM)

runB: $(PROGRAM) cleanB
		@echo "******************* Running ***********************"
		export OMP_NUM_THREADS=$(NUMTHREADS) ; time ./$(PROGRAM)

queue: cleanB
		@echo "*************** Submitting a job ******************"
		echo "cd $$PWD ; export OMP_NUM_THREADS=$(NUMTHREADS) ; time ./$(PROGRAM) " | qsub -q fila -pe threads $(NUMTHREADS) -N job_ldgc
		@echo "***************************************************"

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
		rm -f *.optrpt *.o ldgc *~ fort* *.eco *.con job_ldgc* opt_report.txt time_perf.dat
		@echo "Removing: OK"

cleanB:
		@echo "******************* Cleaning **********************"
		rm -f *~ fort* *.eco *.con job_ldgc* opt_report.txt time_perf.dat
		@echo "Removing: OK"
