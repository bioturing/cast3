INCLUDE_LIB =		-I ../lib/htslib -I ../lib/samtools

CC_FLAG = 		-Wall -pthread -Wfatal-errors -O2 $(INCLUDE_LIB)

HTSLIB =		-L ../lib/htslib

OBJ_MAIN =		bam.o 			\
			readbam.o		\
			utils.o

bam_stats:
	gcc $(CC_FLAG) -c -o utils.o ../lib/utils.c
	gcc $(CC_FLAG) -c -o bam.o ../lib/bam.c
	gcc $(CC_FLAG) -c -o readbam.o readbam.c
	gcc $(CC_FLAG) $(HTSLIB) $(OBJ_MAIN) validate_align.c -o ../bin/validate_align -lz -lhts
	rm -rf *.o
