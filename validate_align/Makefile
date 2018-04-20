CC_FLAG = 		-Wall -pthread -Wfatal-errors -O2

INCLUDE_LIB =		-I ../htslib

HTSLIB =		-L ../htslib

OBJ_MAIN =		readbam.o

bam_stats:
	gcc $(CC_FLAG) $(HTSLIB) $(INCLUDE_LIB) -c -o readbam.o readbam.c
	gcc $(CC_FLAG) $(HTSLIB) $(INCLUDE_LIB) $(OBJ_MAIN) validate_align.c -o ../bin/validate_align -lz -lhts
	rm -rf *.o