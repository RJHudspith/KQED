TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libKQED_cuda.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

LIBFILES=./LIBFILES/all_kernels.c \
	./LIBFILES/cheby.c ./LIBFILES/chnr_dV.c ./LIBFILES/chnr_dS.c \
	./LIBFILES/chnr_dT.c ./LIBFILES/con_kernel.c \
	./LIBFILES/getff-new.c ./LIBFILES/init.c \
	./LIBFILES/io.c ./LIBFILES/kernels.c ./LIBFILES/pi_pert.c \
	./LIBFILES/QED_kernel.c ./LIBFILES/QED_kernel_xy0.c \
	./LIBFILES/SYMXY.c ./LIBFILES/SYMXY0.c \
	./LIBFILES/Tabd.c \
        ./LIBFILES/sub_kernel.c

UTILSFILES=./UTILS/GLU_timer.c ./UTILS/GLU_bswap.c ./UTILS/crc32c.c \
	./UTILS/corr_malloc.c

## all the source files apart from amu_forlattice.c
libKQED_cuda_a_SOURCES = \
	${LIBFILES} \
	${UTILSFILES}
libKQED_cuda_a_AR = $(NVCC) $(NVCCFLAGS) -lib -o
libKQED_cuda_a_LIBADD = $(CUDA_LIBS)
libKQED_cuda_a_CFLAGS = $(CUDA_CFLAGS)

.cu.o:
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $<

## and their respective headers
include_HEADERS = \
	$(LIBFILES:./LIBFILES/%.c=./HEADERS/%.h) \
	$(UTILSFILES:./UTILS/%.c=./HEADERS/%.h) \
	./HEADERS/KQED.h \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = KQED_CUDA

KQED_CUDA_SOURCES = amu_forlattice.c
KQED_CUDA_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
KQED_CUDA_LDADD = libKQED_cuda.a ${LDFLAGS} -lm $(CUDA_LDFLAGS) -lcudart

endif

