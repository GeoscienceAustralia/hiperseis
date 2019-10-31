.PHONY : obj driver clean

#FC = mpifort -mcmodel=medium -shared-intel
#FC = mpifort -mcmodel=medium -O2
#gfortran compilation line:
FC = /opt/mpich/bin/mpifort -O2 -Wall -Wextra -fbounds-check -fbacktrace -pedantic -g -fopenmp

all: run

run: driver

driver: obj RJ_MCMC_RF.f90
	$(FC) -o run RJ_MCMC_RF.f90 \
        theo.o qlayer.o four1.o voro2qmodel.o \
        priorvalue.o whichcell.o dsvdcmp.o \
        svbksb.o dpythag.o -lm

obj: theo.f qlayer.f90 four1.f voro2qmodel.f90 priorvalue.f90 whichcell.f90 dsvdcmp.f svbksb.f dpythag.f
	$(FC) -c theo.f
	$(FC) -c qlayer.f90
	$(FC) -c four1.f
	$(FC) -c voro2qmodel.f90
	$(FC) -c priorvalue.f90
	$(FC) -c whichcell.f90
	$(FC) -c dsvdcmp.f
	$(FC) -c svbksb.f
	$(FC) -c dpythag.f

clean:
	/bin/rm *.o run
