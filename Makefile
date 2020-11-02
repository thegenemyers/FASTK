DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = FastK Histex Tabex Profex Haplex Vennex Homex

all: deflate.lib libhts.a $(ALL)

include HTSLIB/htslib_static.mk

deflate.lib: LIBDEFLATE
	cd LIBDEFLATE; make; cd ..

libhts.a: HTSLIB
	cd HTSLIB; make; cd ..

FastK: FastK.c FastK.h io.c split.c count.c table.c merge.c io.c gene_core.c gene_core.h MSDsort.c LSDsort.c
	gcc $(CFLAGS) -o FastK -I./HTSLIB $(HTSLIB_sstatic_LDFLAGS) FastK.c io.c split.c count.c table.c merge.c MSDsort.c LSDsort.c gene_core.c LIBDEFLATE/libdeflate.a HTSLIB/libhts.a -lpthread $(HTSLIB_static_LIBS)

Histex: Histex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Histex Histex.c gene_core.c -lpthread -lm

Tabex: Tabex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Tabex Tabex.c gene_core.c -lpthread -lm

Profex: Profex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Profex Profex.c gene_core.c -lpthread -lm

Haplex: Haplex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Haplex Haplex.c gene_core.c -lpthread -lm

Vennex: Vennex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Vennex Vennex.c gene_core.c -lpthread -lm

Homex: Homex.c gene_core.c gene_core.h
	gcc $(CFLAGS) -o Homex Homex.c gene_core.c -lpthread -lm

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
