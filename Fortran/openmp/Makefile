# System dependent definitions

FC=	nvfortran
#FFLAGS=	-mp=gpu -Minfo=mp -O3 -gpu=cc80
FFLAGS=	-mp=multicore -Minfo=mp -O3
LFLAGS= -I/sw/spack-levante/nvhpc-22.5-v4oky3/Linux_x86_64/22.5/compilers/include -L/sw/spack-levante/nvhpc-22.5-v4oky3/Linux_x86_64/22.5/compilers/lib -L/sw/spack-levante/nvhpc-22.5-v4oky3/Linux_x86_64/22.5/cuda/11.7/nvvm/lib64 -L/sw/spack-levante/nvhpc-22.5-v4oky3/Linux_x86_64/22.5/cuda/11.7/targets/x86_64-linux/lib -lnvToolsExt

# System independent definitions

MF=	Makefile

EXE=	cfd

SRC= \
        nvtx.f90 \
	boundary.f90 \
	cfd.f90 \
	cfdio.f90 \
	jacobi.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

cfd.o:	boundary.o jacobi.o cfdio.o

tar:
	tar cvf cfd.tar $(MF) $(INC) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) *.mod velocity.dat colourmap.dat cfd.plt core
