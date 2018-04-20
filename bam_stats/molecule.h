#ifndef _MOLECULE_H_
#define _MOLECULE_H_

#include "../lib/utils.h"
#include "../lib/khash.h"

#include "attr.h"

void mlc_init(int n_chr);

void mlc_destroy(int n_chr);

void mlc_insert(int bxid, int pos, int len, struct summary_t *chr_st);

void mlc_get_last(struct summary_t *chr_st);

void gem_merge(struct summary_t *dest, struct summary_t *src);

#endif /* _MOLECULE_H_ */