#include "read_sv.h"

static struct seg_t *retros;
static int n_retro;

static int8_t nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static char nt_char[4] = "ACGT";

/* Read structural variant objects from high confident structural variant file
 * then return number of SVs that have been read
 */
void read_sv(FILE *fp, struct sv_t **svs, int *n_sv, struct genome_t *genome)
{
	size_t len = 0;
	char *s = NULL, alt[1024], ref[1024], name[1024], chr[1024];
	int rc;

	while (getline(&s, &len, fp) != EOF){
		*svs = realloc(*svs, ++(*n_sv) * sizeof(struct sv_t));
		struct sv_t *sv = *svs + *n_sv - 1;
		rc = sscanf(s, "%s\t%d\t%d\t%s\t%s\t%s\t%s",
			    chr, &sv->start, &sv->end, name, ref, alt, sv->type);
		assert(rc == 7 && "Error: Wrong SNP tsv format");

		sv->chr = get_chr_id(chr, genome);
		assert(sv->chr != -1 && "Error: Chr not found in reference");
		sv->cnv = 0;
		sv->seq = NULL;

		/* convert from 1-based coordinate system in VCF file into 0-based 
		 * coordinate system
		 */
		--sv->start, --sv->end;

		/* preprocess some structural variants & get number of duplication */
		if (strcmp(sv->type, "CNV") == 0 || strcmp(sv->type, "DUP") == 0) {
			sv->cnv = alt[3] - '0';
			strcpy(sv->type, "CNV");
		}

		/* just in case deletion of ALU, LINE1, SVA */
		if (strncmp(sv->type, "DEL", 3) == 0) {
			strcpy(sv->type, "CNV");
		} else if (IS_RETRO(sv)) {
			strcpy(sv->type, "INS");
		}

		/* debug read sv */
		// printf("sv %d: %s %d %d %s %s %s %s\n", *n_sv - 1, chr,
		// 	sv->start, sv->end, name, ref, alt, sv->type);
	}
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

void seg2fasta(struct seg_t *seg, int n_seg, int idx, struct genome_t *genome, FILE *fp)
{
	int i, j, k, c;
	fputc('>', fp);
	fputs(genome->ref_name[idx], fp);
	fputc('\n', fp);
	
	for (i = k = 0; i < n_seg; ++i, ++seg) {
		if (seg->strand == '-') {
			for (j = seg->end; j > seg->start; --j) {
				c = seg->seq[j];
				c = nt4_table[c] == 4 ? 4 : 3 - nt4_table[c];
				fputc(nt_char[c], fp);
				if (++k % 80 == 0)
					fputc('\n', fp);
			}
		} else {
			for (j = seg->start; j < seg->end; ++j) {
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
	int sum, i;
	for (i = sum = 0; i < n; ++i, ++seg)
		sum = sum + seg->end - seg->start;
	return sum;
}

/* compare structural variant segment */
static int comp_sv(const void *a, const void *b)
{
	struct sv_t *s1 = (struct sv_t *)a;
	struct sv_t *s2 = (struct sv_t *)b;

	if (s1->chr == s2->chr)
		return s1->start - s2->start;
	return s1->chr - s2->chr;
}

/* This function assumes the input is an sorted structural variant list */
void sim_sv(struct sv_t *svs, int n_sv, struct genome_t *genome, FILE *hap_f, FILE *bed_f)
{
#define PUSH_SEG(_tid, _start, _end, _seq, _strand) do {		\
	segs[_tid] = realloc(segs[_tid], ++seg_cnt[_tid]		\
			     * sizeof(struct seg_t));			\
	segs[_tid][seg_cnt[_tid] - 1] = (struct seg_t){			\
		.start = _start,					\
		.end = _end, 						\
		.seq = _seq,						\
		.strand = _strand					\
	};								\
} while (0);

	int i, tid;
	struct sv_t *sv = svs;
	struct seg_t *segs[genome->sz];
	int seg_cnt[genome->sz], prev[genome->sz];

	for (i = 0; i < genome->sz; ++i) {
		segs[i] = NULL;
		prev[i] = seg_cnt[i] = 0;
	}

	/* sort structural variants by its chromosome then pos */
	qsort(svs, n_sv, sizeof(struct sv_t), comp_sv);

	while (sv - svs < n_sv) {
		tid = sv->chr;

		/* remove overlap regions */
		if (sv->start <= prev[tid]) {
			fprintf(stderr,
				"Warnning: overlap structural variants %s:%d-%d, skipping!\n",
				genome->ref_name[tid], sv->start, sv->end);
			++sv;
			continue;
		}

		PUSH_SEG(tid, prev[tid], sv->start, genome->seq[tid], '+');

		if (IS_DEL(sv)) {
			prev[tid] = sv->end;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tDEL\n", genome->ref_name[tid],
				sv->start, genome->ref_name[tid], sv->end);
		} else if (IS_DUP(sv)) {
			for (i = 0; i < sv->cnv; ++i) {
				PUSH_SEG(tid, sv->start, sv->end, genome->seq[tid], '+');
				fprintf(bed_f, "%s\t%d\t%s\t%d\tDUP\n", genome->ref_name[tid],
					sv->start, genome->ref_name[tid], sv->end);
			}
			prev[tid] = sv->end;
		} else if (IS_INV(sv)) {
			PUSH_SEG(tid, sv->start, sv->end, genome->seq[tid], '-');
			prev[tid] = sv->end;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINV\n", genome->ref_name[tid],
				sv->start, genome->ref_name[tid], sv->end);
		} else if (IS_INS(sv)) {
			double r = ran_normal();
			r *= n_retro;
			int idx = abs((int)(r + 0.5));
			idx = idx > n_retro ? n_retro - 1 : idx;

			PUSH_SEG(tid, retros[idx].start, retros[idx].end,
				 retros[idx].seq, retros[idx].strand);
			prev[tid] = sv->start;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINS\n", genome->ref_name[tid],
				sv->start, genome->ref_name[tid],
				sv->start + (retros[idx].end - retros[idx].start)); 
		} else if (IS_SINS(sv)) {
			PUSH_SEG(tid, sv->start, sv->end, sv->seq, '-');
			prev[tid] = sv->start;
			fprintf(bed_f, "%s\t%d\t%s\t%d\tINS\n", genome->ref_name[tid],
				sv->start, genome->ref_name[tid], sv->end);
		}
		++sv;
	}

	for (i = 0; i < genome->sz; ++i) {
		if (seg_cnt[i] == 0) {
			PUSH_SEG(i, 0, genome->chr_sz[i], genome->seq[i], '+');
		} else {
			PUSH_SEG(i, prev[i], genome->chr_sz[i], genome->seq[i], '+');
		}
		seg2fasta(segs[i], seg_cnt[i], i, genome, hap_f);
	}

#undef PUSH_SEG
}

/* read retrotransposon from static bed file */
void read_retro(FILE *fp, struct genome_t *genome)
{	
	int i = 0, rc, tid;
	size_t len = 0;
	char chr[1024], *s = NULL;

	while (getline(&s, &len, fp) != EOF) {
		retros = realloc(retros, ++n_retro * sizeof(struct seg_t));
		rc = sscanf(s, "%s\t%d\t%d\t%c", chr, &retros[i].start,
			    &retros[i].end, &retros[i].strand);
		assert(rc == 4 && "Error: Wrong retro bed format!");
		tid = get_chr_id(chr, genome);
		assert(tid != -1 && "Error: Chr not found in reference");
		retros[i++].seq = genome->seq[tid];

		/* debug read retro */
		// printf("retro %d: %s, %d, %d, %c\n", i - 1, chr,
		// 	retros[i - 1].start, retros[i - 1].end, retros[i - 1].strand);
	}
}