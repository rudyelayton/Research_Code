TARGET = optical_absorption.x
OBJECTS = optical_absorption.o
MOD_FILES = inputoutput_mod.mod dipolevector_mod.mod readingfile_mod.mod dvec_mod.mod matele_mod.mod
FC = ifort 

FFLAGS =
LDFLAGS =

ifeq (${FC},ifort)
	MKLROOT = /opt/intel/composer_xe_2013_sp1.1.106/mkl
	FFLAGS += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
#	LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmkl_blas95.a
#	LDFLAGS += ${MKLROOT}/lib/ia32/libmkl_lapack95.a
	LDFLAGS += -qopenmp
endif

.SUFFIXES : .o .f90
%.o : %.mod
.f90.o:
	${FC} -c $<

${TARGET} : ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
