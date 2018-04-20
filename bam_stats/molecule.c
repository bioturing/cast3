#include "molecule.h"
#include "khash_barcode.h"

static int *n_mlc;
static struct arr_u64_t **mlc;
static int bx_cnt;
static khash_t(khash_str) *bx_kh;

static void gem_insert(int bxid, int start, int end, int len, int cnt,
		       struct summary_t *chr_st)
{
	if (bxid >= chr_st->n_gem) {
		int old_sz = chr_st->n_gem;
		chr_st->n_gem = bxid + 1;
		chr_st->gem = realloc(chr_st->gem, chr_st->n_gem *
				      sizeof(struct set_mole_t));
		memset(chr_st->gem + old_sz, 0, (chr_st->n_gem - old_sz) *
		       sizeof(struct set_mole_t));
	}

	struct set_mole_t *p = chr_st->gem + bxid;
	p->mlc = realloc(p->mlc, ++p->sz * sizeof(struct mole_t));
	p->mlc[p->sz - 1] = (struct mole_t){
		.start = start, .end = end, .len = len, .cnt = cnt,
		.chr_id = chr_st->chr_id
	};
}

static void do_merge(struct summary_t *dest, int dest_id,
		     struct summary_t *src, int src_id)
{
	if (dest_id == dest->n_gem) {
		dest->gem = realloc(dest->gem, ++dest->n_gem *
				    sizeof(struct set_mole_t));
		dest->gem[dest_id].sz = 0;
		dest->gem[dest_id].mlc = NULL;
		dest->gem[dest_id].barcode = strdup(src->gem[src_id].barcode);
	}

	int new_sz = dest->gem[dest_id].sz + src->gem[src_id].sz;
	dest->gem[dest_id].mlc = realloc(dest->gem[dest_id].mlc, new_sz *
					 sizeof(struct mole_t));
	memcpy(dest->gem[dest_id].mlc + dest->gem[dest_id].sz,
	       src->gem[src_id].mlc, src->gem[src_id].sz * sizeof(struct mole_t));
	dest->gem[dest_id].sz += src->gem[src_id].sz;
}

void gem_merge(struct summary_t *dest, struct summary_t *src)
{
	int i, bxid;

	for (i = 0; i < src->n_gem; ++i) {
		bxid = get_bxid(bx_kh, &bx_cnt, src->barcode[i]);
		src->gem[i].barcode = src->barcode[i];
		do_merge(dest, bxid, src, i);
		free(src->gem[i].mlc);
		free(src->barcode[i]);
	}

	free(src->gem);
	free(src->barcode);
}

void mlc_init(int n_chr)
{
	n_mlc = calloc(n_chr, sizeof(int));
	mlc = calloc(n_chr, sizeof(struct arr*));
	bx_kh = kh_init(khash_str);
}

void mlc_destroy(int n_chr)
{
	int i, j;
	for (i = 0; i < n_chr; ++i) {
		for (j = 0; j < n_mlc[i]; ++j)
			free(mlc[i][j].val);
		free(mlc[i]);
	}
	free(mlc);
	free(n_mlc);
	free_kh_str_bx(bx_kh);
}

void mlc_insert(int bxid, int pos, int len, struct summary_t *chr_st)
{
	int *n = &n_mlc[chr_st->chr_id];
	struct arr_u64_t **p = &mlc[chr_st->chr_id];

	if (bxid >= *n) {
		int old_sz = *n;
		*n = bxid + 1;
		*p = realloc(*p, *n * sizeof(struct arr_u64_t));
		memset(*p + old_sz, 0, (*n - old_sz) * sizeof(struct arr_u64_t));
	}

	struct arr_u64_t *memb = *p + bxid;
	int start, end, mlc_len;

	if (memb->sz &&
	    (memb->val[memb->sz - 1] & MASK32) < pos - MLC_LIMIT_DIS_2READ) {
		end = (memb->val[memb->sz - 1] & MASK32) +
		      (memb->val[memb->sz - 1] >> SHIFT32);
		start = (memb->val[0] & MASK32);
		mlc_len = end - start;
		if (mlc_len >= MIN_MLC_LEN)
			gem_insert(bxid, start, end, mlc_len, memb->sz, chr_st);
		memb->sz = 1;
		memb->val = realloc(memb->val, sizeof(uint64_t));
		memb->val[0] = (1ULL * len << SHIFT32) + pos;
	} else {
		memb->val = realloc(memb->val, ++memb->sz * sizeof(uint64_t));
		memb->val[memb->sz - 1] = (1ULL * len << SHIFT32) + pos;
	}
}

void mlc_get_last(struct summary_t *chr_st)
{
	int n = n_mlc[chr_st->chr_id];
	struct arr_u64_t *p = mlc[chr_st->chr_id];
	int i, start, end, mlc_len;

	for (i = 0; i < n; ++i) {
		struct arr_u64_t *memb = p + i;
		if (memb->sz) {
			end = (memb->val[memb->sz - 1] & MASK32) +
			      (memb->val[memb->sz - 1] >> SHIFT32);
			start = (memb->val[0] & MASK32);
			mlc_len = end - start;
			if (mlc_len >= MIN_MLC_LEN)
				gem_insert(i, start, end, mlc_len, memb->sz, chr_st);
		}
	}
}