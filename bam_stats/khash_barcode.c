#include "khash_barcode.h"

#include <stdio.h>
#include <assert.h>

void free_kh_str_bx(khash_t(khash_str) *bx_kh)
{
	khint_t i;
	for (i = kh_begin(bx_kh); i != kh_end(bx_kh); ++i)
		if (kh_exist(bx_kh, i))
			free((char*)kh_key(bx_kh, i));
	kh_destroy(khash_str, bx_kh);
}

int get_bxid(khash_t(khash_str) *bx_kh, int *bx_khid, char *bx)
{
	int is_new;

	khiter_t k = kh_put(khash_str, bx_kh, bx, &is_new);
	assert(is_new == 0 || is_new == 1);
	if (is_new) {
		kh_value(bx_kh, k) = *bx_khid;
		*bx_khid += 1;
		kh_key(bx_kh, k) = strdup(bx);
	}

	return kh_value(bx_kh, k);
}