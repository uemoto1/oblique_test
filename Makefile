FC=gfortran
LIBS=
INCLUDES=
FFLAGS=-O3 -fopenmp

TARGET=test01

$(TARGET): fdtd1d1.o main.o
	$(FC) -o $@ $^ $(FFLAGS) $(LIBS) 

main.o: main.f90 fdtd1d1.o

fdtd1d1.o: fdtd1d1.f90

%.o: %.f90
	$(FC) -o $@ -c $< $(FFLAGS) $(INCLUDES)




	
