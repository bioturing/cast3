#include "read_stats.h"
#include "molecule.h"
#include "htslib/sam.h"

static int **cell;

void coverage_init(int n_chr, int *chr_len)
{
	int i;
	cell = malloc(n_chr * sizeof(int *));
	for (i = 0; i < n_chr; ++i)
		cell[i] = calloc(chr_len[i] + 1, sizeof(int));
}

void get_read_stats(struct alg_inf_t *read, struct summary_t *chr_st)
{
	int pos, i, isize;

	isize = __abs(read->isize) < N_ISIZE ? __abs(read->isize) : 0;
	++chr_st->isize[isize];

	pos = read->pos;
	for (i = 0; i < read->n_cigar; ++i) {
		int oplen = bam_cigar_oplen(read->cigar[i]);
		char opchr = bam_cigar_opchr(read->cigar[i]);
		if (opchr == 'D') {
			pos += oplen;
		} else if (opchr == 'M') {
			++cell[chr_st->chr_id][pos];
			--cell[chr_st->chr_id][pos + oplen];
			pos += oplen;
		}
	}
}

void coverage_get(struct summary_t *chr_st, int chr_len)
{
	int i, cnt;
	for (i = cnt = 0; i < chr_len; ++i) {
		cnt += cell[chr_st->chr_id][i];
		assert(cnt >= 0);
		if (cnt < N_COVER)
			++chr_st->cover[cnt];
	}
	free(cell[chr_st->chr_id]);
}