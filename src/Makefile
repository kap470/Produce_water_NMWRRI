FC      = gfortran
CPPFLAGS = 
LDFLAGS  = 
FFLAGS  = -g
COMPILE  = $(FC) $(FFLAGS) -c
LINK     = $(FC) -o

OBJECT  = 04_TYPES.o 03_FUNCTIONS.o 02_IO.o 01_Run.o

all: program

program: $(OBJECT)
	$(LINK) program.exe $(OBJECT)

%.o: ./%.f90
	$(COMPILE) $<

clean-all: clean
	rm -f *.exe

clean:
	rm -f *.o *.mod *.mod0
