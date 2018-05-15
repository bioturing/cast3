#include "htslib/khash.h"

#include "attr.h"
#include "duplicate.h"
#include "molecule.h"
#include "read_stats.h"
#include "khash_barcode.h"

struct prog_args {
	int nthread;
	int distance_thres;
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
	return -((*(struct pair_t *)a).first - (*(struct pair_t *)b).first);
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
	char path[BUFSZ];
	FILE *fi_sum, *fi_cplot, *fi_mplot, *fi_mlc;

	sprintf(path, "%s/summary.inf", args.outdir);
	if (!(fi_sum = fopen(path, "w")))
		__PERROR("Could not open summary.inf");
	sprintf(path, "%s/cover_plot.tsv", args.outdir);
	if (!(fi_cplot = fopen(path, "w")))
		__PERROR("Could not open cover_plot.tsv");
	sprintf(path, "%s/mlc_plot.tsv", args.outdir);
	if (!(fi_mplot = fopen(path, "w")))
		__PERROR("Could not open mlc_plot.tsv");
	sprintf(path, "%s/molecule.tsv", args.outdir);
	if (!(fi_mlc = fopen(path, "w")))
		__PERROR("Could not open molecule.tsv");

	int i, j, total_gem_detected = 0;
	int64_t total_dna_20kb = 0, total_dna_100kb = 0;
	int64_t total_mlc_len = 0, total_mlc_cnt = 0, n_mlc = 0;
	struct arr_pair_t mlc = (struct arr_pair_t){.val = NULL, .sz = 0};
	int64_t mlc_plot[N_MLC];
	memset(mlc_plot, 0, N_MLC * sizeof(int64_t));

	fprintf(fi_mlc, "#chr\tstart\tend\tbarcode\tread_count\n");

	for (i = 0; i < genome_st.n_gem; ++i) {
		if (!genome_st.gem[i].sz) {
			free(genome_st.gem[i].barcode);
			continue;
		}
		++total_gem_detected;
		n_mlc += genome_st.gem[i].sz;

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
			mlc.val[mlc.sz - 1] = (struct pair_t){.first = cnt, .second = len};

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
		int cnt = mlc.val[i].first;
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
	fprintf(fi_sum, "Mean Molecule per GEM:\t%ld\n",
		n_mlc / total_gem_detected);
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
	samFile *bam_f;
	if (!(bam_f = sam_open(file_name, "r")))
		__PERROR("Could not open BAM file");
	if (!(bam_inf.bam_i = sam_index_load(bam_f, file_name)))
		__ERROR("BAM file must be indexed!\n");

	/* Init the header */
	bam_inf.b_hdr = sam_hdr_read(bam_f);

	sam_close(bam_f);
	return bam_inf;
}

static void read_bam_unmapped(struct bam_inf_t *bam_inf)
{
	samFile *bam_f = sam_open(bam_inf->bam_path, "r");
	hts_itr_t *iter = sam_itr_queryi(bam_inf->bam_i, HTS_IDX_NOCOOR, 0, 0);
	bam1_t *b = bam_init1();
	int i;

	/* all bam record here is unmapped */
	while (sam_itr_next(bam_f, iter, b) >= 0) {
		++genome_st.total_read;
		++genome_st.total_unmap;
		genome_st.total_len += b->core.l_qseq;

		uint8_t *qual = bam_get_qual(b);
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

static void read_bam_target(hts_itr_t *iter, samFile *bam_f,
			    struct summary_t *chr_st)
{
	bam1_t *b = bam_init1();
	uint8_t *tag_data;
	char *bar_s;
	int bx_khid = 0, n_alg_inf = 0, bxid, i;
	struct alg_inf_t *alg_inf = NULL;
	khash_t(khash_str) *bx_kh = kh_init(khash_str);

	while (sam_itr_next(bam_f, iter, b) >= 0) {
		/* skip align is supplementary or not primary */
		if (b->core.flag & (FLAG_NOT_PRI | FLAG_SUPPLEMENT))
			continue;	
		++chr_st->total_read;
		chr_st->total_len += b->core.l_qseq;

		uint8_t *qual = bam_get_qual(b);
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

		if (b->core.n_cigar == 0) {
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
	samFile *bam_f = sam_open(bam_inf->bam_path, "r");
	hts_itr_t *iter;
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
			__ERROR("Genome file is not correct!\n");
		struct summary_t chr_st;
		memset(&chr_st, 0, sizeof(struct summary_t));
		chr_st.chr_id = id;
		iter = sam_itr_queryi(bam_inf->bam_i, id, 1,
				      bam_inf->b_hdr->target_len[id]);
		read_bam_target(iter, bam_f, &chr_st);
		__VERBOSE("Done chr %s in %.2fs.\n",
			  bam_inf->b_hdr->target_name[id],
			  realtime() - local_time);
	} while (1);

	sam_close(bam_f);
}

static void read_bam(struct bam_inf_t *bam_inf)
{
	old_time = realtime();
	__VERBOSE("Get bam stats ... \n");

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

	__VERBOSE("Done all in %.2fs.\n", realtime() - old_time);

	/* output result */
	__VERBOSE("Output result ... \n");
	output_result(bam_inf);

	/* free memory allocate */
	mlc_destroy(bam_inf->b_hdr->n_targets);
	free(pthr);
	pthread_attr_destroy(&attr);
	pthread_mutex_destroy(&lock_id);
	pthread_mutex_destroy(&lock_merge);

	__VERBOSE("Finish!\n");
}

static void load_genome(char *file_path)
{
	__VERBOSE("Loading genome ... \n");
	old_time = realtime();

	int i, id, c;
	size_t len = 0;
	char *s = NULL;
	FILE *fi;
	if (!(fi = fopen(file_path, "r")))
		__PERROR("Could not open genome file!\n");

	c = getc(fi);
	if (c != '>')
		__ERROR("Genome file is not fasta format!\n");
	ungetc(c, fi);

	chr_len = NULL;
	id = -1;
	while (getline(&s, &len, fi) != EOF) {
		// __VERBOSE("%s", s);
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
	__VERBOSE("Done in %.2fs.\n", realtime() - old_time);
}

static void print_usage()
{
	__VERBOSE("./bstats <option> genome_file bam_file\n");
	__VERBOSE("\n");
	__VERBOSE("Version: %s\n", VERSION);
	__VERBOSE("\n");
	__VERBOSE("Input bam_file must be indexed and must has BX tag.\n");
	__VERBOSE("Molecules were detected by cluster read into group.\n");
	__VERBOSE("In each group, distance between two consecutive reads lower than -d.\n");
	__VERBOSE("Filter molecule has length < 1000 or number of reads on it < 4.\n");
	__VERBOSE("\n");
	__VERBOSE("Option:\n");
	__VERBOSE("  -t INT                number of threads [1]\n");
	__VERBOSE("  -o STR                output directory [./]\n");
	__VERBOSE("  -d STR                threshold for distance between two consecutive\n");
	__VERBOSE("                        reads in a molecule [50000]\n");
	__VERBOSE("\n");
	__VERBOSE("This tool will generate some file:\n");
	__VERBOSE("  summary.inf           summary of bam file\n");
	__VERBOSE("  mlc_plot.tsv          molecule length histogram\n");
	__VERBOSE("  cover_plot.tsv        coverage histogram\n");
	__VERBOSE("  molecule.tsv          all molecule info\n");
}

static void get_args(int argc, char *argv[])
{
	int c;
	args.outdir = "./";
	args.nthread = 1;
	args.distance_thres = DIS_THRES_DEFAULT;

	while ((c = getopt(argc, argv, "sht:o:d:")) >= 0) {
		switch (c) {
		case 't':
			args.nthread = atoi(optarg);
			break;
		case 'o':
			args.outdir = optarg;
			break;
		case 'd':
			args.distance_thres = atoi(optarg);
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
	make_dir(args.outdir);
}

int main(int argc, char *argv[])
{
	get_args(argc, argv);
	load_genome(args.genome_path);
	struct bam_inf_t bam_inf = init_bam_reader(args.bam_path);
	mlc_init(bam_inf.b_hdr->n_targets, args.distance_thres);
	coverage_init(bam_inf.b_hdr->n_targets, chr_len);
	read_bam(&bam_inf);

	return 0;
}