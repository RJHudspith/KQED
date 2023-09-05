TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libKQED_cuda.a

## just include the header directory
# AM_CFLAGS=-I${TOPDIR}/src_cuda/HEADERS/
AM_CPPFLAGS=-I./HEADERS/

LIBFILES=./LIBFILES/all_kernels.cu \
	./LIBFILES/cheby.cu ./LIBFILES/chnr_dV.cu ./LIBFILES/chnr_dS.cu \
	./LIBFILES/chnr_dT.cu ./LIBFILES/con_kernel.cu \
	./LIBFILES/getff-new.cu ./LIBFILES/init.cu \
	./LIBFILES/io.cu ./LIBFILES/kernels.cu ./LIBFILES/pi_pert.cu \
	./LIBFILES/QED_kernel.cu ./LIBFILES/QED_kernel_xy0.cu \
	./LIBFILES/SYMXY.cu ./LIBFILES/SYMXY0.cu \
	./LIBFILES/Tabd.cu \
        ./LIBFILES/sub_kernel.cu

UTILSFILES=./UTILS/GLU_timer.c ./UTILS/GLU_bswap.c ./UTILS/crc32c.c

## all the source files apart from amu_forlattice.c
libKQED_cuda_a_SOURCES = ${LIBFILES} ${UTILSFILES}
libKQED_cuda_a_AR = $(NVCC) $(NVCCFLAGS) -lib -o
libKQED_cuda_a_LIBADD = $(CUDA_LIBS)
libKQED_cuda_a_CFLAGS = $(CUDA_CFLAGS) $(INCLUDES) $(CFLAGS)

.cu.o:
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) $(DEFS) $(DEFAULT_INCLUDES) $(CFLAGS) $(CPPFLAGS) $(AM_CPPFLAGS) -dc -o $@ $<
.c.o:
	$(NVCC) -x cu $(NVCCFLAGS) $(INCLUDES) $(DEFS) $(DEFAULT_INCLUDES) $(CFLAGS) $(CPPFLAGS) $(AM_CPPFLAGS) -dc -o $@ $<

## and their respective headers
include_HEADERS = \
	$(LIBFILES:./LIBFILES/%.c=./HEADERS/%.h) \
	$(UTILSFILES:./UTILS/%.c=./HEADERS/%.h) \
	./HEADERS/KQED.h \
	./HEADERS/cu_util.h \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = KQED_CUDA

KQED_CUDA_SOURCES = amu_forlattice.cu
KQED_CUDA_CFLAGS = $(CUDA_CFLAGS) $(INCLUDES) $(CFLAGS)
KQED_CUDA_LDADD = libKQED_cuda.a ${LDFLAGS} -lm $(CUDA_LIBS) $(CUDA_LDFLAGS) -lcudart
KQED_CUDA_LINK = $(NVCC) $(NVCCFLAGS) $(KQED_CUDA_CFLAGS) $(KQED_CUDA_LDFLAGS) $(LDFLAGS) -o $@

endif

