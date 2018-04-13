#ifndef _READ_FASTA_H_
#define _READ_FASTA_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

struct genome_t {
	char **seq;
	char **ref_name;
	int *chr_sz;
	int sz;
};

struct seq_t {
	char *name;
	char *seq;
};

int get_chr_id(char *s, struct genome_t *genome);

struct genome_t read_genome(FILE *fp);

#endif /* _READ_FASTA_H_ */

