## everything is in src or src_cuda
if USE_CUDA
SUBDIRS = src_cuda
else
SUBDIRS = src
endif

ACLOCAL_AMFLAGS = -I m4

TOPDIR=@top_srcdir@
REG_FILES = \
	$(TOPDIR)/PRECOMP/REGexample2.txt \
	$(TOPDIR)/PRECOMP/REGexample4.txt \
	$(TOPDIR)/PRECOMP/REGexample5.txt \
	$(TOPDIR)/PRECOMP/REGexample6.txt \
	$(TOPDIR)/PRECOMP/REGexample7.txt \
	$(TOPDIR)/PRECOMP/REGexample10.txt \
	$(TOPDIR)/PRECOMP/REGress.sh
bin_SCRIPTS = $(REG_FILES)

PRECOMP_FILES = \
	$(TOPDIR)/PRECOMP/taylorx_cksum.bin \
	$(TOPDIR)/PRECOMP/taylory_cksum.bin \
	$(TOPDIR)/PRECOMP/FFxy_single_cksum.bin
precompdir = $(prefix)/PRECOMP
precomp_DATA = $(PRECOMP_FILES)

