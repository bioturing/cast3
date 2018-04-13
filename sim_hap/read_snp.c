#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "read_snp.h"
#include "read_sv.h"

#define MAX_SNP_LINE			2000
#define MAX_INDEL_SEQ_BUFFER		1000

void read_snp(FILE *fp, struct sv_t *svs, int *n, char **genome, FILE *bed_f)
{
	int i;
	size_t l;
	char *line = malloc(MAX_SNP_LINE);	
	char *ref = malloc(MAX_INDEL_SEQ_BUFFER);
	char *alt = malloc(MAX_INDEL_SEQ_BUFFER);
	int start, end;
	char chr[3];
	struct sv_t *sv = svs + *n;
	extern char *chroms[];

	while (getline(&line, &l, fp) != EOF) {
		sscanf(line, "%s%d%d%s%s", chr, &start, &end, ref, alt);
		/* convert to 0-based coordinate */
		--start, --end;
		CHROM_IDX(chr, i);
		
		if (end - start + 1 > MAX_INDEL_SEQ_BUFFER) {
			fprintf(stderr, "Too long short-INDEL!!!\n");
			exit(EXIT_FAILURE);
		}

		if (start - end == 0) {
			genome[i][start] = alt[0];
		} else {
			++*n;
			//printf("Indel %d\n", j++);
			assert(end > start);
			*(sv) = (struct sv_t){.start = start, .end = end, .chr = i};
			strcpy(sv->type, strlen(ref) > strlen(alt) ? "CNV" : "SINS");
			if (!strcmp(sv->type, "INS")) {
				sv->seq = malloc(strlen(alt));
				strcpy(sv->seq, alt);
			}
			++sv;
		}
	}
}
