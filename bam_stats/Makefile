INCLUDE_LIB =		-I ../lib/htslib -I ../lib/samtools

CC_FLAG = 		-Wall -pthread -Wfatal-errors -O2 $(INCLUDE_LIB)

HTSLIB =		-L ../lib/htslib

OBJ_MAIN =		bam.o 			\
			molecule.o		\
			duplicate.o		\
			read_stats.o		\
			khash_barcode.o		\
			utils.o

bam_stats:
	gcc $(CC_FLAG) -c -o utils.o ../lib/utils.c
	gcc $(CC_FLAG) -c -o bam.o ../lib/bam.c
	gcc $(CC_FLAG) -c -o molecule.o molecule.c
	gcc $(CC_FLAG) -c -o khash_barcode.o khash_barcode.c
	gcc $(CC_FLAG) -c -o duplicate.o duplicate.c
	gcc $(CC_FLAG) -c -o read_stats.o read_stats.c
	gcc $(CC_FLAG) $(HTSLIB) $(OBJ_MAIN) bam_stats.c -o ../bin/bstats -lz -lhts
	rm -rf *.o
