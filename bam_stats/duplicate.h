#ifndef _DUPLICATE_H_
#define _DUPLICATE_H_

#include "attr.h"

void checkdup_process(struct summary_t *chr_st, int n_alg_inf,
		      struct alg_inf_t *alg_inf);

void checkdup_insert(bam1_t *b, int bxid, struct summary_t *chr_st, int *n_alg_inf,
		     struct alg_inf_t **alg_inf);

#endif /* _DUPLICATE_H_ */