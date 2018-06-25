FC=gfortran
LIBS=
INCLUDES=
FFLAGS=-O3
TARGET=test01

$(TARGET): inputoutput.o fdtd1d1.o main.o 
	$(FC) $^ -o $@ $(FFLAGS) $(LIBS) 

inputoutput.f90: inputoutput.py
	python inputoutput.py

main.o: main.f90 fdtd1d1.mod

inputoutput.o: inputoutput.f90
inputoutput.mod: inputoutput.f90 inputoutput.o

fdtd1d1.o: fdtd1d1.f90 inputoutput.mod
fdtd1d1.mod: fdtd1d1.f90 fdtd1d1.o

	



.PHONY: clean all

all: $(TARGET)
	
clean:
	rm *.o $(TARGET)


%.o: %.f90
	$(FC) $< -c -o $@ $(FFLAGS)

%.mod: %.f90 %.o
	@true


	
