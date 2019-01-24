FC = gfortran
FLAG = -O2
#FC           = /opt/ibmcmp/xlf/8.1/bin/xlf
#FLAG         = -O3 -qstrict -qfree=f90 -qxlf77=leadzero -qarch=auto -qmaxmem=16384 -qtune=auto -qsuffix=f=f90
#FLAG         =  -CA -CB -CS -CU -CV -WB
executable   = polar_tensors_traj.exe
 
ALL: main

OBJECTS = main.o kset.o dcell.o  erfunc.o  invert.o \
          readin.o realE_dippim.o recipE_dippim.o \
          setup.o  conjgradpim.o  modules.o boxreset.o \
          polarizability.o wannier2dipoles.o	

main : $(OBJECTS)
	$(FC) -o $(executable) $(OBJECTS) $(FLAG)

clean :
	rm -f *.o *.mod

.SUFFIXES : .f90
.f90.o:
	$(FC) -c $(FLAG) $*.f90

# dependencies

main.o : main.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod readin.o setup.o \
         conjgradpim.o polarizability.o wannier2dipoles.o

kset.o : kset.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod dcell.o

dcell.o : dcell.f90 BOXDATA.mod

erfunc.o : erfunc.f90 

readin.o : readin.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod dcell.o 

realE_dippim.o : realE_dippim.f90 COMMONDATA.mod BOXDATA.mod

recipE_dippim.o : recipE_dippim.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod \
                  erfunc.o

setup.o : setup.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod dcell.o invert.o boxreset.o kset.o 

conjgradpim.o : conjgradpim.f90 COMMONDATA.mod BOXDATA.mod RECIPDATA.mod \
                erfunc.o 

boxreset.o : boxreset.f90 COMMONDATA.mod BOXDATA.mod dcell.o invert.o

modules.o : modules.f90

polarizability.o : polarizability.f90 invert.o COMMONDATA.mod

wannier2dipoles.o : wannier2dipoles.f90 COMMONDATA.mod BOXDATA.mod

COMMONDATA.mod : modules.f90
	$(FC) -c $(FLAG) modules.f90

BOXDATA.mod : modules.f90
	$(FC) -c $(FLAG) modules.f90

RECIPDATA.mod : modules.f90
	$(FC) -c $(FLAG) modules.f90
