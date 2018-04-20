#include "../lib/khash.h"

#include "attr.h"
#include "duplicate.h"
#include "molecule.h"
#include "read_stats.h"
#include "khash_barcode.h"

/* how to get align element:
 * name:     %s, bam1_qname(b)
 *
 * flag:     %d, b->core.flag
 *
 * tid:      %d, b->core.tid (its name get from bam header)
 *
 * pos:      %d, b->core.pos
 *
 * qual:     %d, b->core.qual
 *
 * cigar:    uint32_t *cigar = bam1_cigar(b) 
 *           %d%c, bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i])
 *
 * mtid:     %d, b->core.mtid (its name get from bam header)
 *
 * mpos:     %d, b->core.mpos
 *
 * ins_sz:   %d, b->core.isize
 *
 * seq:      uint8_t *seq = bam1_seq(b)
 *           %c, bam_nt16_rev_table[bam1_seqi(seq, i)]
 *
 * qual_str: uint8_t *qual = bam1_qual(b)
 *           %c, qual[i] + 33
 *
 * barcode:  tag_data = bam_aux_get(b, "BX");
 *           bar_s = bam_aux2Z(tag_data);
 */

struct prog_args {
	int nthread;
	char *bam_path;
	char *genome_path;
	char *outdir;
} args;

static double old_time;
static int *chr_len;
static int64_t sum_nt4, sum_nt_amb;
static struct summary_t genome_st;

static pthread_mutex_t lock_merge = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t lock_id = PTHREAD_MUTEX_INITIALIZER;

static inline int cmpfunc_int(const void *a, const void *b)
{
	return *(int *)a - *(int *)b;
}

static inline int cmpfunc_mlc(const void *a, const void *b)
{
	return -((*(uint64_t *)a >> SHIFT32) - (*(uint64_t *)b >> SHIFT32));
}

static void merge_to_genome(struct summary_t *chr_st)
{
	genome_st.total_read += chr_st->total_read;
	genome_st.total_dup += chr_st->total_dup;
	genome_st.total_unmap += chr_st->total_unmap;
	genome_st.total_len += chr_st->total_len;
	genome_st.q30base_r1 += chr_st->q30base_r1;
	genome_st.q30base_r2 += chr_st->q30base_r2;
	genome_st.nbase_r1 += chr_st->nbase_r1;
	genome_st.nbase_r2 += chr_st->nbase_r2;

	int i;
	for (i = 0; i < N_COVER; ++i)
		genome_st.cover[i] += chr_st->cover[i];
	for (i = 0; i < N_ISIZE; ++i)
		genome_st.isize[i] += chr_st->isize[i];

	gem_merge(&genome_st, chr_st);
}

static int get_median_isize(struct summary_t *p)
{
	int i;
	int64_t sum = 0;
	for (i = 1; i < N_ISIZE; ++i) {
		sum += p->isize[i];
		if (sum >= (p->total_read >> 1))
			return i;
	}
	return N_ISIZE;
}

static void output_result(struct bam_inf_t *bam_inf)
{
	genome_st.cover[0] -= sum_nt_amb;
	char path[MAX_DIR_LEN];
	FILE *fi_sum, *fi_cplot, *fi_mplot, *fi_mlc;

	sprintf(path, "%s/summary.inf", args.outdir);
	if (!(fi_sum = fopen(path, "w")))
		__perror("Could not open summary.inf");
	sprintf(path, "%s/cover_plot.tsv", args.outdir);
	if (!(fi_cplot = fopen(path, "w")))
		__perror("Could not open cover_plot.tsv");
	sprintf(path, "%s/mlc_plot.tsv", args.outdir);
	if (!(fi_mplot = fopen(path, "w")))
		__perror("Could not open mlc_plot.tsv");
	sprintf(path, "%s/molecule.tsv", args.outdir);
	if (!(fi_mlc = fopen(path, "w")))
		__perror("Could not open molecule.tsv");

	int i, j, total_gem_detected = 0;
	int64_t total_dna_20kb = 0, total_dna_100kb = 0;
	int64_t total_mlc_len = 0, total_mlc_cnt = 0;
	struct arr_u64_t mlc = (struct arr_u64_t){.val = NULL, .sz = 0};
	int64_t mlc_plot[N_MLC];
	memset(mlc_plot, 0, N_MLC * sizeof(int64_t));

	fprintf(fi_mlc, "#chr\tstart\tend\tbarcode\tread_count\n");

	for (i = 0; i < genome_st.n_gem; ++i) {
		if (!genome_st.gem[i].sz) {
			free(genome_st.gem[i].barcode);
			continue;
		}
		++total_gem_detected;

		for (j = 0; j < genome_st.gem[i].sz; ++j) {
			struct mole_t *cur_mlc = genome_st.gem[i].mlc + j;
			int len = cur_mlc->len;
			int cnt = cur_mlc->cnt;
			int start = cur_mlc->start;
			int end = cur_mlc->end;
			int chr_id = cur_mlc->chr_id;
			char *barcode = genome_st.gem[i].barcode;
			fprintf(fi_mlc, "%s\t%d\t%d\t%s\t%d\n",
				bam_inf->b_hdr->target_name[chr_id],
				start, end, barcode, cnt);
			if (len >= MLC_20KB)
				total_dna_20kb += len;
			if (len >= MLC_100KB)
				total_dna_100kb += len;
			total_mlc_len += len;
			total_mlc_cnt += cnt;
			++mlc.sz;
			if ((mlc.sz & (mlc.sz - 1)) == 0)
				mlc.val = realloc(mlc.val, (mlc.sz << 1) *
						  sizeof(uint64_t));
			mlc.val[mlc.sz - 1] = (1ULL * cnt << SHIFT32) + len;

			int bin_id = len / MLC_BIN_PLOT;
			if (bin_id < N_MLC)
				mlc_plot[bin_id] += len;
		}

		free(genome_st.gem[i].barcode);
	}

	qsort(mlc.val, mlc.sz, sizeof(uint64_t), cmpfunc_mlc);
	int64_t sum = 0;
	int n50_read_per_mlc = -1;
	for (i = sum = 0; i < mlc.sz; ++i) {
		int cnt = (int)(mlc.val[i] >> SHIFT32);
		sum += cnt;
		if (sum >= (total_mlc_cnt >> 1)) {
			n50_read_per_mlc = cnt;
			break;
		}
	}
	assert(n50_read_per_mlc != -1);

	/* output summary data */
	if (genome_st.nbase_r1 != 0)
		fprintf(fi_sum, "Mode: Paired-end\n");
	else
		fprintf(fi_sum, "Mode: Single-end\n");
	fprintf(fi_sum, "Number of Reads:\t%lu\n", genome_st.total_read);
	fprintf(fi_sum, "PCR Duplication:\t%.2f%%\n",
		1.0 * genome_st.total_dup / genome_st.total_read * 100);
	fprintf(fi_sum, "Mapped Reads:\t%.2f%%\n",
		1.0 * (genome_st.total_read - genome_st.total_unmap) /
		genome_st.total_read * 100);
	if (genome_st.nbase_r1 != 0) {
		fprintf(fi_sum, "Q30 bases, Read 1:\t%.2f%%\n",
			1.0 * genome_st.q30base_r1 / genome_st.nbase_r1 * 100);
		fprintf(fi_sum, "Q30 bases, Read 2:\t%.2f%%\n",
			1.0 * genome_st.q30base_r2 / genome_st.nbase_r2 * 100);
	} else {
		fprintf(fi_sum, "Q30 bases:\t%.2f%%\n",
			1.0 * genome_st.q30base_r2 / genome_st.nbase_r2 * 100);
	}
	fprintf(fi_sum, "Zero Coverage:\t%.2f%%\n",
		1.0 * genome_st.cover[0] / sum_nt4 * 100);
	if (genome_st.nbase_r1 != 0)
		fprintf(fi_sum, "Median Insert Size:\t%d\n", get_median_isize(&genome_st));
	fprintf(fi_sum, "Mean Depth:\t%.1fX\n", 1.0 * genome_st.total_len / sum_nt4);
	fprintf(fi_sum, "GEMs Detected:\t%d\n", total_gem_detected);
	fprintf(fi_sum, "Mean DNA per GEM:\t%ld\n",
		total_mlc_len / total_gem_detected);
	fprintf(fi_sum, "DNA in Molecules >20kb:\t%.1f%%\n",
		1.0 * total_dna_20kb / total_mlc_len * 100);
	fprintf(fi_sum, "DNA in Molecules >100kb:\t%.1f%%\n",
		1.0 * total_dna_100kb / total_mlc_len * 100);
	fprintf(fi_sum, "Mean Molecule Length:\t%ld\n",
		total_mlc_len / mlc.sz);
	fprintf(fi_sum, "N50 Reads per Molecule (LPM):\t%d\n", n50_read_per_mlc);

	/* output data into 2 plot */
	for (i = 0; i < N_COVER; ++i)
		fprintf(fi_cplot, "%d\t%lu\n", i, genome_st.cover[i]);
	for (i = 0; i < N_MLC; ++i)
		fprintf(fi_mplot, "%d\t%lu\n", i * MLC_BIN_PLOT / 1000, mlc_plot[i]);

	fclose(fi_sum);
	fclose(fi_cplot);
	fclose(fi_mplot);
	fclose(fi_mlc);
}

struct bam_inf_t init_bam_reader(char *file_name)
{
	struct bam_inf_t bam_inf;
	bam_inf.bam_path = file_name;
	bam_inf.cur_id = 0;

	/* Open bam file and the index file */
	bamFile bam_f;
	if (!(bam_f = bam_open(file_name, "r")))
		__perror("Could not open BAM file");
	if (!(bam_inf.bam_i = bam_index_load(file_name)))
		__error("BAM file must be indexed!");

	/* Init the header */
	bam_inf.b_hdr = bam_header_read(bam_f);

	bam_close(bam_f);
	return bam_inf;
}

static void read_bam_unmapped(struct bam_inf_t *bam_inf)
{
	bamFile bam_f = bam_open(bam_inf->bam_path, "r");
	bam_iter_t iter = bam_iter_query(bam_inf->bam_i, HTS_IDX_NOCOOR, 0, 0);
	bam1_t *b = bam_init1();
	int i;

	/* all bam record here is unmapped */
	while (bam_iter_read(bam_f, iter, b) >= 0) {
		++genome_st.total_read;
		++genome_st.total_unmap;
		genome_st.total_len += b->core.l_qseq;

		uint8_t *qual = bam1_qual(b);
		assert((b->core.flag & (FLAG_NOT_PRI | FLAG_SUPPLEMENT)) == 0);
		if (b->core.flag & FLAG_READ1) {
			for (i = 0; i < b->core.l_qseq; ++i) {
				++genome_st.nbase_r1;
				if (qual[i] >= 30)
					++genome_st.q30base_r1;
			}
		} else {
			for (i = 0; i < b->core.l_qseq; ++i) {
				++genome_st.nbase_r2;
				if (qual[i] >= 30)
					++genome_st.q30base_r2;
			}
		}	
	}
}

static void read_bam_target(bam_iter_t iter, bamFile bam_f,
			    struct summary_t *chr_st)
{
	bam1_t *b = bam_init1();
	uint8_t *tag_data;
	char *bar_s;
	int bx_khid = 0, n_alg_inf = 0, bxid, i;
	struct alg_inf_t *alg_inf = NULL;
	khash_t(khash_str) *bx_kh = kh_init(khash_str);

	while (bam_iter_read(bam_f, iter, b) >= 0) {
		/* skip align is supplementary or not primary */
		if (b->core.flag & (FLAG_NOT_PRI | FLAG_SUPPLEMENT))
			continue;	
		++chr_st->total_read;
		chr_st->total_len += b->core.l_qseq;

		uint8_t *qual = bam1_qual(b);
		if (b->core.flag & FLAG_READ1) {
	 		for (i = 0; i < b->core.l_qseq; ++i) {
	 			++chr_st->nbase_r1;
	 			if (qual[i] >= 30)
	 				++chr_st->q30base_r1;
	 		}
	 	} else {
	 		for (i = 0; i < b->core.l_qseq; ++i) {
	 			++chr_st->nbase_r2;
	 			if (qual[i] >= 30)
	 				++chr_st->q30base_r2;
	 		}
	 	}

		if (!b->core.n_cigar) {
			++chr_st->total_unmap;
			continue;
		}

		tag_data = bam_aux_get(b, "BX");
		assert(tag_data);
		bar_s = bam_aux2Z(tag_data);
		bxid = get_bxid(bx_kh, &bx_khid, bar_s);
		checkdup_insert(b, bxid, chr_st, &n_alg_inf, &alg_inf);

		assert(bxid <= chr_st->n_barcode);
		if (bxid == chr_st->n_barcode) {
			chr_st->barcode = realloc(chr_st->barcode,
						  ++chr_st->n_barcode *
						  sizeof(char *));
			chr_st->barcode[bxid] = strdup(bar_s);
		}
	}

	checkdup_process(chr_st, n_alg_inf, alg_inf);
	coverage_get(chr_st, chr_len[chr_st->chr_id]);
	mlc_get_last(chr_st);

	pthread_mutex_lock(&lock_merge);
	merge_to_genome(chr_st);
	pthread_mutex_unlock(&lock_merge);

	free_kh_str_bx(bx_kh);
	free(alg_inf);
}

static void *process(void *data)
{
	struct bam_inf_t *bam_inf = (struct bam_inf_t *)data;
	bamFile bam_f = bam_open(bam_inf->bam_path, "r");
	bam_iter_t iter;
	int id;

	do {
		pthread_mutex_lock(&lock_id);
		if (bam_inf->cur_id == bam_inf->b_hdr->n_targets) {
			pthread_mutex_unlock(&lock_id);
			pthread_exit(NULL);
			return NULL;
		} else {
			id = bam_inf->cur_id++;
		}
		pthread_mutex_unlock(&lock_id);

		double local_time = realtime();
		if (bam_inf->b_hdr->target_len[id] != chr_len[id])
			__ferror("Genome file is not correct!");
		struct summary_t chr_st;
		memset(&chr_st, 0, sizeof(struct summary_t));
		chr_st.chr_id = id;
		iter = bam_iter_query(bam_inf->bam_i, id, 1,
				      bam_inf->b_hdr->target_len[id]);
		read_bam_target(iter, bam_f, &chr_st);
		__verbose("Done chr %s in %.2fs.\n",
			  bam_inf->b_hdr->target_name[id],
			  realtime() - local_time);
	} while (1);

	bam_close(bam_f);
}

static void read_bam(struct bam_inf_t *bam_inf)
{
	old_time = realtime();
	__verbose("Get bam stats ... \n");

	/* init data */
	mlc_init(bam_inf->b_hdr->n_targets);
	coverage_init(bam_inf->b_hdr->n_targets, chr_len);

	/* Get unmapped read, tid == -1 */
	read_bam_unmapped(bam_inf);

	/* get bam stats by parallel */
	pthread_t *pthr = calloc(args.nthread, sizeof(pthread_t));;
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	int i;
	for (i = 0; i < args.nthread; ++i)
		pthread_create(&pthr[i], &attr, process, bam_inf);
	for (i = 0; i < args.nthread; ++i)
		pthread_join(pthr[i], NULL);

	__verbose("Done all in %.2fs.\n", realtime() - old_time);

	/* output result */
	__verbose("Output result ... \n");
	output_result(bam_inf);

	/* free memory allocate */
	mlc_destroy(bam_inf->b_hdr->n_targets);
	free(pthr);
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&lock_id);
	pthread_mutex_destroy(&lock_merge);

	__verbose("Finish!\n");
}

static void load_genome(char *file_path)
{
	__verbose("Loading genome ... \n");
	old_time = realtime();

	int i, id, c;
	size_t len = 0;
	char *s = NULL;
	FILE *fi;
	if (!(fi = fopen(file_path, "r")))
		__perror("Could not open genome file");

	c = getc(fi);
	if (c != '>')
		__error("Genome file is not fasta format!\n");
	ungetc(c, fi);

	chr_len = NULL;
	id = -1;
	while (getline(&s, &len, fi) != EOF) {
		// __verbose("%s", s);
		if (s[0] == '>') {
			++id;
			chr_len = realloc(chr_len, (id + 1) * sizeof(int));
			chr_len[id] = 0;
			continue;
		}
		for (i = 0; i < len && s[i] != '\n'; ++i) {
			++chr_len[id];
			if (s[i] != 'N' && s[i] != 'n')
				++sum_nt4;
		}
	}

	for (i = 0; i <= id; ++i)
		sum_nt_amb += chr_len[i];
	sum_nt_amb -= sum_nt4;

	fclose(fi);
	__verbose("Done in %.2fs.\n", realtime() - old_time);
}

static void print_usage()
{
	__verbose("./bstats <option> genome_file bam_file\n");
	__verbose("\n");
	__verbose("input bam_file is sorted by position and must be index\n");
	__verbose("\n");
	__verbose("option:\n");
	__verbose("  -t INT           number of threads [1]\n");
	__verbose("  -o STR           output directory [./]\n");
	__verbose("\n");
	__verbose("this tool will generate some file:\n");
	__verbose("  summary.inf      summary of bam file\n");
	__verbose("  mlc_plot.tsv          molecule length histogram\n");
	__verbose("  cover_plot.tsv        coverage histogram\n");
	__verbose("  molecule.tsv     all molecule info\n");
}

static void get_args(int argc, char *argv[])
{
	int c;
	args.outdir = "./";
	args.nthread = 1;

	while ((c = getopt(argc, argv, "sht:o:")) >= 0) {
		switch (c) {
		case 't':
			args.nthread = atoi(optarg);
			break;
		case 'o':
			args.outdir = optarg;
			break; 
		case 'h':
			print_usage();
			exit(0);
		default:
			print_usage();
			exit(1);
		}
	}

	if (optind + 2 == argc) {
		args.genome_path = argv[optind];
		args.bam_path = argv[optind + 1];
	} else {
		print_usage();
		exit(1);
	}

	args.nthread = __min(args.nthread, 16);
	make_outdir(args.outdir);
}

int main(int argc, char *argv[])
{
	get_args(argc, argv);
	load_genome(args.genome_path);
	struct bam_inf_t bam_inf = init_bam_reader(args.bam_path);
	read_bam(&bam_inf);

	return 0;
}