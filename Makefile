DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = FastK Fastrm Fastmv Fastcp Histex Tabex Profex Haplex Homex Vennex Logex

all: deflate.lib libhts.a $(ALL)

include HTSLIB/htslib_static.mk

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

libhts.a: HTSLIB
	cd HTSLIB; make; cd ..

libfastk.c : gene_core.c
libfastk.h : gene_core.h

FastK: FastK.c FastK.h io.c split.c count.c table.c merge.c io.c gene_core.c gene_core.h MSDsort.c LSDsort.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o FastK -I./HTSLIB $(HTSLIB_static_LDFLAGS) FastK.c io.c split.c count.c table.c merge.c MSDsort.c LSDsort.c libfastk.c LIBDEFLATE/libdeflate.a HTSLIB/libhts.a -lpthread $(HTSLIB_static_LIBS)

Fastrm: Fastrm.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Fastrm Fastrm.c gene_core.c -lpthread -lm

Fastmv: Fastxfer.c gene_core.c gene_core.h
	gcc $(CFLAGS) -DMOVE -o Fastmv Fastxfer.c gene_core.c -lpthread -lm

Fastcp: Fastxfer.c gene_core.c gene_core.h
	gcc $(CFLAGS) -UMOVE -o Fastcp Fastxfer.c gene_core.c -lpthread -lm

Histex: Histex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Histex Histex.c libfastk.c -lpthread -lm

Tabex: Tabex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Tabex Tabex.c libfastk.c -lpthread -lm

Profex: Profex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Profex Profex.c libfastk.c -lpthread -lm

Haplex: Haplex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Haplex Haplex.c libfastk.c -lpthread -lm

Homex: Homex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Homex Homex.c libfastk.c -lpthread -lm

Vennex: Vennex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Vennex Vennex.c libfastk.c -lpthread -lm

Logex: Logex.c libfastk.c libfastk.h
	gcc $(CFLAGS) -o Logex Logex.c libfastk.c -lpthread -lm

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
