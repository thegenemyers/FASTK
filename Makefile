DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FastK Fastrm Fastmv Fastcp Fastmerge Histex Tabex Profex Logex Vennex Symmex Haplex Homex Fastcat

all: deflate.lib libhts.a $(ALL)

include HTSLIB/htslib_static.mk

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

libhts.a: HTSLIB
	cd HTSLIB; make libhts.a; cd ..

HTSLIB/htslib_static.mk:
	cd HTSLIB; make htslib_static.mk; cd ..

libfastk.c : gene_core.c
libfastk.h : gene_core.h

FastK: FastK.c FastK.h io.c split.c count.c table.c merge.c io.c gene_core.c gene_core.h MSDsort.c LSDsort.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o FastK -I./HTSLIB $(HTSLIB_static_LDFLAGS) FastK.c io.c split.c count.c table.c merge.c MSDsort.c LSDsort.c libfastk.c LIBDEFLATE/libdeflate.a HTSLIB/libhts.a -lpthread $(HTSLIB_static_LIBS)

Fastrm: Fastrm.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o Fastrm Fastrm.c gene_core.c -lpthread -lm

Fastmv: Fastxfer.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -DMOVE -o Fastmv Fastxfer.c gene_core.c -lpthread -lm

Fastcp: Fastxfer.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -UMOVE -o Fastcp Fastxfer.c gene_core.c -lpthread -lm

Fastmerge: Fastmerge.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Fastmerge Fastmerge.c libfastk.c -lpthread -lm

Fastcat: Fastcat.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Fastcat Fastcat.c libfastk.c -lpthread -lm

Histex: Histex.c libfastk.c libfastk.h ONElib.h ONElib.c
	$(CC) $(CFLAGS) -o Histex Histex.c libfastk.c ONElib.c -lpthread -lm

Tabex: Tabex.c libfastk.c libfastk.h ONElib.h ONElib.c
	$(CC) $(CFLAGS) -o Tabex Tabex.c libfastk.c ONElib.c -lpthread -lm

Profex: Profex.c libfastk.c libfastk.h ONElib.h ONElib.c
	$(CC) $(CFLAGS) -o Profex Profex.c libfastk.c ONElib.c -lpthread -lm

Logex: Logex.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Logex Logex.c libfastk.c -lpthread -lm

Vennex: Vennex.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Vennex Vennex.c libfastk.c -lpthread -lm

Symmex: Symmex.c libfastk.c libfastk.h LSDsort.c
	$(CC) $(CFLAGS) -o Symmex Symmex.c libfastk.c LSDsort.c -lpthread -lm

Haplex: Haplex.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Haplex Haplex.c libfastk.c -lpthread -lm

Homex: Homex.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Homex Homex.c libfastk.c -lpthread -lm


tidyup:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastK.tar.gz

clean:
	cd LIBDEFLATE; make clean; cd ..
	cd HTSLIB; make clean; cd ..
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastK.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastK.tar.gz LICENSE README.md Makefile *.h *.c
