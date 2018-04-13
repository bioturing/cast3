#include "read_fasta.h"

#define SEQ_BUFFER 			20000000
#define NAME_BUFFER 			10
#define MAX_CHROMS 			200
#define MAX_NAME_BUFFER 		200

#define CHROM_IDX(chr, i) \
	for (i = 0; i < MAX_CHROMS && strcmp(chr, *(chroms + i)); i++);

char *chroms[MAX_CHROMS];
int nchroms;

const uint8_t nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void get_name(char *s, FILE *stream)
{
	int c, i = 0;
	while (!isspace(c = getc(stream))) {
		s[i++] = c;
	}
	s[i] = '\0';
	while (c != '\n')
		c = getc(stream);
}

int get_seq(char **s, FILE *stream)
{
	int c, n = 0, buff_i = 1;
	if ((*s = malloc(buff_i*SEQ_BUFFER)) == NULL) {
		fprintf(stderr, "Can not allocate memmory!\n");
		exit(EXIT_FAILURE);
	}
	c = getc(stream);
	while (c != '>' && c != EOF) {
		if (n > SEQ_BUFFER * buff_i) {
			*s = realloc(*s, ++buff_i * SEQ_BUFFER);
		}
		if (c != '\n')
			(*s)[n++] = c;
		c = getc(stream);
	}
	ungetc(c, stream);
	return n;
}

char **read_genome(FILE *fp)
{
	int i = 0, c, n;
	char name[NAME_BUFFER];
	char **seqs = (char **)calloc(MAX_CHROMS, sizeof(char *));
	extern char *chroms[];
	extern int nchroms;

	for (i = 0; i < MAX_CHROMS; i++){
		chroms[i] = malloc(MAX_NAME_BUFFER);
	}

	i = 0;
	while ((c = getc(fp)) != EOF){
		if (c == '>'){
			get_name(name, fp);
			i = i >= MAX_CHROMS?MAX_CHROMS - 1:i;
			strcpy(chroms[i], name);
			n = get_seq(seqs + i, fp);
			++i;
		}
	}
	nchroms= i;
	return(seqs);
}