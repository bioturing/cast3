CC =			gcc

CFLAGS =		-Wall -Wfatal-errors -O2

OBJS = 			read_fasta.o 			\
			read_snp.o 			\
			read_sv.o

MAIN_OBJ =		../bin/sim.exe

sim_hap:
	$(CC) $(CFLAGS) -c -o read_fasta.o read_fasta.c
	$(CC) $(CFLAGS) -c -o read_snp.o read_snp.c
	$(CC) $(CFLAGS) -c -o read_sv.o read_sv.c
	$(CC) $(CFLAGS) $(OBJS) main.c -o $(MAIN_OBJ) -lm
	rm -rf *.o
