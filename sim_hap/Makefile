CC =			gcc

CFLAGS =		-static -Wall -Wfatal-errors -O2

LINK =			-lm

OBJS = 			read_fasta.o 			\
			read_snp.o

MAIN_OBJ =		../bin/sim

all:
	$(CC) $(CFLAGS) -c -o read_fasta.o read_fasta.c $(LINK)
	$(CC) $(CFLAGS) -c -o read_snp.o read_snp.c $(LINK)
	$(CC) $(CFLAGS) $(OBJS) read_sv.c -o $(MAIN_OBJ) $(LINK)
clean:
	rm read_sv
