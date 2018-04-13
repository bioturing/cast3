#ifndef READ_FASTA_H
#define READ_FASTA_H

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#define MAX_CHROMS 200

void get_name(char *s, FILE *stream);
int get_seq(char **s, FILE *stream);
char **read_genome(FILE *fp);
void read_fasta(FILE *fp, char **seqs);

char *chroms[MAX_CHROMS];
int nchroms;
const uint8_t nst_nt4_table[256];

struct seq_t{
	char *name;
	char *seq;
};

#endif

