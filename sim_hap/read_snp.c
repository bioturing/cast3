#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "read_snp.h"

#define LIMIT_SNP_LEN			1000

void read_snp(FILE *fp, struct sv_t **svs, int *n_sv, struct genome_t *genome)
{
	int tid, rc, start, end;
	size_t len = 0;
	char *line = NULL;
	char ref[1024], alt[1024], chr[1024];

	while (getline(&line, &len, fp) != EOF) {
		rc = sscanf(line, "%s\t%d\t%d\t%s\t%s", chr, &start, &end, ref, alt);
		assert(rc == 5 && "Error: Wrong SNP tsv format");
		assert(end - start + 1 < LIMIT_SNP_LEN &&
		       "Error: Too long short-INDEL!!!");

		/* convert to 0-based coordinate */
		--start, --end;
		tid = get_chr_id(chr, genome);
		assert(tid != -1 && "Error: Chr not found in reference");

		if (end - start == 0) {
			assert(start < genome->chr_sz[tid] &&
			       "Error: SNP overflow");
			genome->seq[tid][start] = alt[0];
		} else {
			assert(end > start && "Error: SV has end <= start");
			++(*n_sv);
			*svs = realloc(*svs, *n_sv * sizeof(struct sv_t));
			struct sv_t *sv = *svs + *n_sv - 1;
			*sv = (struct sv_t){.start = start, .end = end, .chr = tid};
			strcpy(sv->type, strlen(ref) > strlen(alt) ? "CNV" : "SINS");
			sv->cnv = 0;
			if (strcmp(sv->type, "SINS") == 0)
				sv->seq = strdup(alt);
		}
	}
}
