#ifndef _ATTR_H_
#define _ATTR_H_

#include <stdint.h>

#include "htslib/sam.h"

/* Sam flag */
#define FLAG_PAIR		0x001	/* Read has its mate */
#define FLAG_PROPER 		0x002	/* Read and its mate are proper pair */
#define FLAG_UNMAP 		0x004	/* Read is unmapped */
#define FLAG_M_UNMAP		0x008	/* Read's mate is unmapped */
#define FLAG_REVERSE 		0x010	/* Read is reversed */
#define FLAG_M_REVERSE 		0x020	/* Read's mate is reversed */
#define FLAG_READ1 		0x040	/* Read is first read */
#define FLAG_READ2 		0x080	/* Read is second read */
#define FLAG_NOT_PRI		0x100	/* Alginment is not primary */
#define FLAG_DUPLICATE 		0x400	/* Read is pcr or duplicated */
#define FLAG_SUPPLEMENT		0x800	/* Alignment is supplementary */

/* summary attribute */
#define N_COVER			200
#define N_MLC			100
#define N_ISIZE 		1000

#define MLC_LIMIT_DIS_2READ	100000
#define MIN_MLC_LEN		4000
#define MLC_20KB		20000
#define MLC_100KB		100000
#define MLC_BIN_PLOT		4000

struct bam_inf_t {
	char *bam_path;
	hts_idx_t *bam_i;
	bam_hdr_t *b_hdr;
	int cur_id;
};

struct mole_t {
	int start;
	int end;
	int len;
	int cnt;
	int chr_id;
};

struct set_mole_t {
	struct mole_t *mlc;
	char *barcode;
	int sz;
};

struct arr_u64_t {
	uint64_t *val;
	int sz;
};

struct summary_t {
	int chr_id;
	int64_t total_read;
	int64_t total_dup;
	int64_t total_unmap;
	int64_t total_len;
	int64_t q30base_r1;
	int64_t q30base_r2;
	int64_t nbase_r1;
	int64_t nbase_r2;
	int64_t isize[N_ISIZE];
	int64_t cover[N_COVER];
	struct set_mole_t *gem;
	int n_gem;
	char **barcode;
	int n_barcode;
};

struct alg_inf_t {
	int mtid;
	int mpos;
	int bxid;
	int tid;
	int pos;
	int flag;
	int isize;
	int len;
	uint32_t *cigar;
	int n_cigar;
};

#endif /* _ATTR_H_ */