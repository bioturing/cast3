#include "read_fasta.h"

#define SEQ_BUFFER 			1048576

/* return -1 if not found */
int get_chr_id(char *s, struct genome_t *genome)
{
	int i;
	for (i = 0; i < genome->sz; ++i)
		if (strcmp(genome->ref_name[i], s) == 0)
			return i;
	return -1;
}

/* get chr name upto space or endline */
static char *get_name(FILE *stream)
{
	char *s = malloc(1);
	int sz = 1, i = 0, c;
	while (1) {
		c = getc(stream);
		if (isspace(c)) // isspace(): \t, \n, \v, \f, \r
			break;
		if (i + 1 == sz) {
			sz <<= 1;
			s = realloc(s, sz);
		}
		s[i++] = c;
	}
	s = realloc(s, i + 1);
	s[i] = '\0';
	while (c != '\n')
		c = getc(stream);
	return s;
}

static char *get_seq(FILE *stream, int *sz)
{
	char *s = malloc(SEQ_BUFFER);
	char *line;
	int buf_cnt = 1;
	int c;
	size_t l_buff = 0, ls = 0, l = 0;

	while ((c = fgetc(stream)) != EOF){
		/* roll back 1 character */
		ungetc(c, stream);
		if (c == '>')
			break;
		l = getline(&line, &l_buff, stream);
		if (ls + l > buf_cnt * SEQ_BUFFER){
			buf_cnt <<= 1;
			s = realloc(s, buf_cnt * SEQ_BUFFER);
		}
		memccpy(s + ls, line, '\n', l);
		ls += l - 1; //exclude the '\n'
		s[ls] = '\0';
	}
	*sz = ls;
	return s;
}


struct genome_t read_genome(FILE *fp)
{
	struct genome_t genome;
	genome.sz = 0;
	genome.seq = genome.ref_name = NULL;
	genome.chr_sz = NULL;

	while (1) {
		int c = getc(fp);
		if (c == EOF)
			break;
		if (c == '>') {
			++genome.sz;
			genome.seq = realloc(genome.seq, genome.sz * sizeof(char *));
			genome.ref_name = realloc(genome.ref_name, genome.sz * sizeof(char *));
			genome.chr_sz = realloc(genome.chr_sz, genome.sz * sizeof(int));
			genome.ref_name[genome.sz - 1] = get_name(fp);
			genome.seq[genome.sz - 1] = get_seq(fp, &genome.chr_sz[genome.sz - 1]);
		} else {
			fprintf(stderr, "Genome file is not fasta format!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	return genome;
}
