#include "readbam.h"

static samFile *b_f;
static bam_hdr_t *b_hdr;

int get_tid(char *s)
{
	assert(s);
	return bam_name2id(b_hdr, s);
}

char *get_name(int tid)
{
	if (tid == -1)
		return "*";
	return b_hdr->target_name[tid];
}

void init_bam(char *file_path)
{
	if (!(b_f = hts_open(file_path, "r"))) {
		fprintf(stderr, "Error: Cannot open BAM file: %s", file_path);
		exit(EXIT_FAILURE);
	}
	b_hdr = sam_hdr_read(b_f);
}

char *convert_scigar(uint32_t *bcigar, int sz)
{
	char *ret = NULL, tmp[1024];
	int i, len = 0, tlen;

	for (i = 0; i < sz; ++i) {
		sprintf(tmp, "%d%c", bam_cigar_oplen(bcigar[i]),
			bam_cigar_opchr(bcigar[i]));
		tlen = strlen(tmp);
		ret = realloc(ret, (len + tlen) * sizeof(char));
		memcpy(ret + len, tmp, tlen);
		len += tlen;
	}
	ret = realloc(ret, (len + 1) * sizeof(char));
	ret[len] = '\0';

	return ret;
}

void sam_write(FILE *fp, bam1_t *b)
{
	uint8_t *tag_data = bam_aux_get(b, "BX");
	if (!tag_data) {
		fprintf(stderr, "Error: Please add barcode info to your sam record (tag BX:Z)\n");
		exit(EXIT_FAILURE);
	}
	char *bar_code = bam_aux2Z(tag_data);
	char *alternative = NULL;
	tag_data = bam_aux_get(b, "XA");
	if (tag_data)
		alternative = bam_aux2Z(tag_data);
	fprintf(fp, "%d\t%s\t%d\t%s\t%s\t%s\n",
		b->core.flag, get_name(b->core.tid), b->core.pos + 1,
		convert_scigar(bam_get_cigar(b), b->core.n_cigar),
		bar_code, alternative);
}

int read_bam(bam1_t *b)
{
	while (1) {
		int rc = sam_read1(b_f, b_hdr, b);
		if (rc <= 0)
			return 0;
		if (b->core.flag & ((1 << 11) + (1 << 8)))
			continue;
		break;
	}
	return 1;
}