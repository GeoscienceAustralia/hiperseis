.PHONY : obj driver clean

# You must make sure mpifort is available on path, e.g. by 'module load openmpi/1.10.2'
ifeq ($(TARGET_ARCH),AVX2)
  # This build is much faster on Intel CPUs with AVX2 instruction set.
  # To generate this build on NCI, run 'setup_make_avx2.sh' shell script first.
  FC = mpifort -mcmodel=medium -O2 -march=core-avx2 -axCORE-AVX2 -xCORE-AVX2
  # Uncomment to add runtime checks option (slows down runtime)
#  FC += -check all
  TARGET = run_avx2
else
  FC = mpifort -mcmodel=medium -O2
  # Uncomment to add runtime checks option (slows down runtime)
#  FC += -fcheck-array-temporaries -fbounds-check
  TARGET = run
endif

all: run

run: driver

driver: obj RJ_MCMC_RF.f90
	$(FC) -o $(TARGET) RJ_MCMC_RF.f90 \
        theo.o qlayer.o four1.o voro2qmodel.o \
        priorvalue.o whichcell.o dsvdcmp.o \
        svbksb.o dpythag.o -lm

obj: theo.f qlayer.f four1.f voro2qmodel.f90 priorvalue.f90 whichcell.f90 dsvdcmp.f svbksb.f dpythag.f
	$(FC) -c theo.f
	$(FC) -c qlayer.f
	$(FC) -c four1.f
	$(FC) -c voro2qmodel.f90
	$(FC) -c priorvalue.f90
	$(FC) -c whichcell.f90
	$(FC) -c dsvdcmp.f
	$(FC) -c svbksb.f
	$(FC) -c dpythag.f

clean:
	/bin/rm *.o run
