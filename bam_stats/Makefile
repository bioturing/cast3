CC_FLAG = 		-Wall -pthread -Wfatal-errors -O2

INCLUDE_LIB =		-I ../htslib

HTSLIB =		-L ../htslib

OBJ_MAIN =		molecule.o		\
			duplicate.o		\
			read_stats.o		\
			khash_barcode.o		\
			utils.o

bam_stats:
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) -c -o utils.o utils.c
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) -c -o molecule.o molecule.c
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) -c -o khash_barcode.o khash_barcode.c
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) -c -o duplicate.o duplicate.c
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) -c -o read_stats.o read_stats.c
	gcc $(CC_FLAG) $(INCLUDE_LIB) $(HTSLIB) $(OBJ_MAIN) bam_stats.c -o ../bin/bstats -lz -lhts
	rm -rf *.o
