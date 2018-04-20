#ifndef _KHASH_BARCODE_H_
#define _KHASH_BARCODE_H_

#include "../lib/khash.h"

KHASH_MAP_INIT_STR(khash_str, int)

void free_kh_str_bx(khash_t(khash_str) *bx_kh);

int get_bxid(khash_t(khash_str) *bx_kh, int *bx_khid, char *bx);

#endif /* _KHASH_BARCODE_H_ */