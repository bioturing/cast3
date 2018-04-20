#include "duplicate.h"
#include "read_stats.h"
#include "molecule.h"

static inline int cmpfunc_alg(const void *a, const void *b)
{
	return ((struct alg_inf_t *)a)->bxid -
	       ((struct alg_inf_t *)b)->bxid;
}

static void assign_alg(bam1_t *b, int bxid, struct alg_inf_t *p)
{
	p->mtid = b->core.mtid;
	p->mpos = b->core.mpos;
	p->tid = b->core.tid;
	p->pos = b->core.pos;
	p->mtid = b->core.mtid;
	p->flag = b->core.flag;
	p->isize = b->core.isize;
	p->n_cigar = b->core.n_cigar;
	p->len = b->core.l_qseq;
	p->bxid = bxid;
	p->cigar = malloc(b->core.n_cigar * sizeof(uint32_t));
	memcpy(p->cigar, bam_get_cigar(b), b->core.n_cigar * sizeof(uint32_t));
}

void checkdup_process(struct summary_t *chr_st, int n_alg_inf,
		      struct alg_inf_t *alg_inf)
{
	int u = 0, v = 0;
	int i, j;

	qsort(alg_inf, n_alg_inf, sizeof(struct alg_inf_t), cmpfunc_alg);

	while (u < n_alg_inf) {
		get_read_stats(&alg_inf[u], chr_st);
		mlc_insert(alg_inf[u].bxid, alg_inf[u].pos, alg_inf[u].len, chr_st);
		while (v < n_alg_inf && alg_inf[v].bxid == alg_inf[u].bxid)
			++v;

		for (i = u + 1; i < v; ++i) {
			for (j = u; j < i; ++j)
				if (alg_inf[i].mpos == alg_inf[j].mpos &&
				    alg_inf[i].mtid == alg_inf[j].mtid)
					break;
			if (j < i) {
				++chr_st->total_dup;
			} else {
				get_read_stats(&alg_inf[i], chr_st);
				mlc_insert(alg_inf[i].bxid, alg_inf[i].pos,
					   alg_inf[i].len, chr_st);
			}
		}
		u = v;
	}

	for (i = 0; i < n_alg_inf; ++i)
		free(alg_inf[i].cigar);
}

void checkdup_insert(bam1_t *b, int bxid, struct summary_t *chr_st,
		     int *n_alg_inf, struct alg_inf_t **alg_inf)
{
	if (*n_alg_inf && b->core.pos != (*alg_inf)[0].pos) {
		checkdup_process(chr_st, *n_alg_inf, *alg_inf);
		*n_alg_inf = 1;
		*alg_inf = realloc(*alg_inf, sizeof(struct alg_inf_t));
		assign_alg(b, bxid, &(*alg_inf)[0]);
		return;
	}

	*alg_inf = realloc(*alg_inf, ++(*n_alg_inf) * sizeof(struct alg_inf_t));
	assign_alg(b, bxid, &(*alg_inf)[*n_alg_inf - 1]);
}