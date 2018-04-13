#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "read_sv.h"

static struct seg_t retros[MAX_RETROS];
static int NUMBER_RETROS = 0;

/* Read structural variant objects from high confident structural variant file
 * then return number of SVs that have been read
*/
void read_sv(FILE *fp, struct sv_t *svs, int *n)
{
	size_t l;
	char *s = malloc(2*BUFFER_SIZE);
	char alt[10], ref[BUFFER_SIZE], name[NAME_BUFFER], chr[3];
	struct sv_t *sv = svs;
	extern char *chroms[];

	while (getline(&s, &l, fp) != EOF){
		sscanf(s, "%s%d%d%s%s%s%s",
		       chr, &sv->start, &sv->end, name, ref, alt, sv->type);
		CHROM_IDX(chr, sv->chr);
		// convert from 1-based coordinate system in VCF file into 0-based 
		// coordinate system
		--sv->start;
		--sv->end;
		// preprocess some structural variants
		// get number of duplication
		if (!strcmp(sv->type, "CNV") || !strcmp(sv->type, "DUP")) {
			sv->cnv = GET_CNV(alt);
			strcpy(sv->type, "CNV");
		}
		// just in case deletion of ALU, LINE1, SVA
		if (!strncmp(sv->type, "DEL", 3)) {
			strcpy(sv->type, "CNV");
		} 
		else if (IS_RETRO(sv)) {
			strcpy(sv->type, "INS");
		}
	
		// Keep the first base for deletion
		/* if (IS_DEL(sv))
			++sv->start;
		*/
		//printf("type: %s, cnv: %d\n", sv->type, sv->cnv);
		sv++;	
	}
	*n = sv-svs;
}

double ran_normal()
{ 
	static int iset = 0; 
	static double gset; 
	double fac, rsq, v1, v2; 
	if (iset == 0) {
		do { 
			v1 = 2.0 * drand48() - 1.0;
			v2 = 2.0 * drand48() - 1.0; 
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq); 
		gset = v1 * fac; 
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

void seg2fasta(struct seg_t *seg, int idx, int n, char ** genome, FILE *fp)
{
	int i, j, k = 0;
	int c, c1;
	char nt_char[4] = "ACGT";
	fputc('>', fp);
	fputs(chroms[idx], fp);
	fputc('\n', fp);
	
	for (i = 0; i < n; i++, seg++) {
		if (seg->strand == '-') {
			for (j = seg->end; j > seg->start; j--) {
			    c = seg->seq[j];
				c1 = nst_nt4_table[c] > 3 ? 4: 3 - nst_nt4_table[c];
				fputc(nt_char[c1], fp);
				if (++k % 80 == 0)
					fputc('\n', fp);
			}
		} else {
			for (j = seg->start; j < seg->end; j++) {
				fputc(seg->seq[j], fp);
				if (++k % 80 == 0)
					fputc('\n', fp);
			}

		}
	}

	if (k % 80 != 0)
		fputc('\n', fp);
}

/* compute sum size of segments */
int segsize(struct seg_t *seg, int n)
{
	int s = 0;
	int i;
	for (i = 0; i < n; i++, seg++) {
		s = s + (seg->end - seg->start);
	}
	return s;
}
/* compare structural variant segment */
int comp_sv(const void *a, const void *b)
{
	struct sv_t *s1 = (struct sv_t *)a;
	struct sv_t *s2 = (struct sv_t *)b;

	if (s1->chr == s2->chr)
		return s1->start - s2->start;
	return s1->chr - s2->chr;
}

void sim_sv(struct sv_t *svs, int n, char **seq, struct seg_t **segs, int *n_segs,
	    FILE *stream, FILE *bed_f)
{
	 /* This function assumes the input is an sorted structural variant list */
	//printf("Number of structural variants: %d\n", n);
	extern int nchroms;
	int i, ch;
	struct sv_t *sv = svs;
	struct seg_t *seg_prev[nchroms];
	int prev[nchroms];
	double r;
	int idx;

	for (i = 0 ; i < nchroms; i++) {
		segs[i] = calloc(MAX_DUP * n, sizeof(struct seg_t));
		prev[i] = 0;
	}
	memmove(seg_prev, segs, nchroms * sizeof(struct seg_t *));
	/* sort structural variants by its pos and chromosome */
	qsort(svs, n, sizeof(struct sv_t), comp_sv);

	while (sv - svs < n) {
		//make one segment chain for each chromosome
		ch = sv->chr;
		//avoid overlap regions
		if (sv->start <= prev[ch]) {
			fprintf(stderr, "WARNNING: overlaid structural variants %s:%d-%d\n",
					chroms[ch], sv->start, sv->end);
			++sv;
			continue;
			//exit(EXIT_FAILURE);
		}
		*(segs[ch]++) = (struct seg_t){.start = prev[ch], .end = sv->start, 
				.seq = *(seq + sv->chr), .strand = '+'};
		//printf("%d\t%d\tNOR\n", prev[ch], sv->start);

		if (IS_DEL(sv)) {
			prev[ch] = sv->end;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tDEL\n", chroms[ch],
				sv->start, chroms[ch], sv->end);
		} else if (IS_DUP(sv)) {
			//printf("cnv: %d\n", sv->cnv);
			for (i = 0; i < sv->cnv; i++) {
				*(segs[ch]++) = (struct seg_t){.start = sv->start, .end = sv->end, 
						.seq = *(seq + sv->chr), .strand = '+'};
				fprintf(bed_f, "%s\t%d\t%s\t%d\tDUP\n",
					chroms[ch], sv->start, chroms[ch], sv->end);
			}
			prev[ch] = sv->end;
		} else if (IS_INV(sv)) {
			*(segs[ch]++) = (struct seg_t){.start = sv->start, .end = sv->end, 
					.seq = *(seq + sv->chr), .strand = '-'};
			prev[ch] = sv->end;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINV\n", chroms[ch],
				sv->start, chroms[ch], sv->end);
		} else if (IS_INS(sv)) {
			r = ran_normal();
			r *= NUMBER_RETROS;
			idx = abs((int)(r + 0.5));
			idx = idx > NUMBER_RETROS?NUMBER_RETROS - 1 : idx;

			*(segs[ch]++) = (struct seg_t){
				.start = retros[idx].start,
				.end = retros[idx].end,
				.seq = retros[idx].seq,
				.strand = retros[idx].strand
			};
			prev[ch] = sv->start;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINS\n", chroms[ch], sv->start, 
				chroms[ch], sv->start + (retros[idx].end - retros[idx].start)); 
		} else if (IS_SINS(sv)) {
			*(segs[ch]++) = (struct seg_t) {
				.start = sv->start,
				.end = sv->end,
				.seq = sv->seq,
				.strand = '+'
			};
			prev[ch] = sv->start;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINS\n", chroms[ch],
				sv->start, chroms[ch], sv->end);
		}
		++sv;
	}
	//Only generate sequences for chromosomes 1-22 X
	for (i = 0 ; i < nchroms - 1; i++) {
		// for chromosome doesn't have any SV
		if (segs[i] == seg_prev[i]) {
			*(segs[i]) = (struct seg_t){.start = 0, 
				.end = strlen(*(seq + i)), .seq = *(seq + i)};
		} else {
			*(segs[i]) = (struct seg_t){.start = (segs[i]-1)->end, 
				.end = strlen(*(seq + i)), .seq = *(seq + i)};
		}
	//	printf("Total size of chrom %s: %d\n", chroms[i], segsize(seg_prev[i], 
	//					segs[i] - seg_prev[i] + 1));
		seg2fasta(seg_prev[i], i, segs[i] - seg_prev[i] + 1, seq, stream);
	}
}

/* read retrotransposon from static bed file */

void read_retro(FILE *fp, char **seqs)
{	
	int c, i = 0, chr_i = 0;
	size_t l;
	char chr[3];
	char *s = malloc(NAME_BUFFER);
	extern char *chroms[];

	while ((c = getline(&s, &l, fp)) != EOF) {
		sscanf(s, "%s\t%d\t%d\t%c", chr, &retros[i].start, &retros[i].end,
			&retros[i].strand);
		CHROM_IDX(chr, chr_i);
		retros[i].seq = *(seqs + chr_i);
		//printf("retro %d:, %s, %d, %d, %c\n", i, chr, retros[i].start, retros[i].end,
		//		retros[i].strand);
		++i;
		++NUMBER_RETROS;
	}
}

int main(int argc, char *argv[])
{
	FILE *sv_fp = fopen(argv[1], "r");
	FILE *genome_fp = fopen(argv[2], "r");
	FILE *retro_fp = fopen(argv[3], "r");
	FILE *snp_fp = fopen(argv[4], "r");
	FILE *bed_f = fopen(argv[5], "w");
	FILE *stream = stdout;
	struct sv_t *svs = (struct sv_t*)calloc(MAX_SV, sizeof(struct sv_t));
	extern int nchroms;
	if (sv_fp == NULL)
		exit(EXIT_FAILURE);
	
	char **genome = read_genome(genome_fp);

	/* initiate some variables related to number of chromosome */
	struct seg_t *segs[nchroms]; 
	int n_sv;
	int n_segs[nchroms];

	read_sv(sv_fp, svs, &n_sv);
	//printf("Done reading genome!\n");
	read_retro(retro_fp, genome);

	read_snp(snp_fp, svs, &n_sv, genome, bed_f);
	//printf("Done read sv!\n");
	sim_sv(svs, n_sv, genome, segs, n_segs, stream, bed_f);
	return 0;
}