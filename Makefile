FC = gfortran

#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 
FFLAGS = -g

DRIVER1_90 = qpexampledriver.f90

LIBRARY = qp.f90
LAPACK = -llapack
BLAS = -lblas

all :  qpspecial

qpspecial : $(DRIVER1_90) $(LIBRARY) $(LAPACK) $(BLAS)
	$(FC) -o qpspecial $(FFLAGS) $(DRIVER1_90) $(LIBRARY) -L/usr/local/lib $(LAPACK) $(BLAS)
