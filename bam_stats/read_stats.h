#ifndef _READ_STATS_H_
#define _READ_STATS_H_

#include "../lib/utils.h"

#include "attr.h"

void coverage_init(int n_chr, int *chr_len);

void get_read_stats(struct alg_inf_t *read, struct summary_t *chr_st);

void coverage_get(struct summary_t *chr_st, int chr_len);

#endif /* _READ_STATS_H_ */