#include "read_fasta.h"
#include "read_sv.h"
#include "read_snp.h"

int main(int argc, char *argv[])
{
	FILE *sv_fp = fopen(argv[1], "r");
	if (!sv_fp) {
		fprintf(stderr, "Could not open sv file!\n");
		exit(EXIT_FAILURE);
	}
	FILE *genome_fp = fopen(argv[2], "r");
	if (!genome_fp) {
		fprintf(stderr, "Could not open genome file!\n");
		exit(EXIT_FAILURE);
	}
	FILE *retro_fp = fopen(argv[3], "r");
	if (!retro_fp) {
		fprintf(stderr, "Could not open retro file!\n");
		exit(EXIT_FAILURE);
	}
	FILE *snp_fp = fopen(argv[4], "r");
	if (!snp_fp) {
		fprintf(stderr, "Could not open snp file!\n");
		exit(EXIT_FAILURE);
	}
	char bed[1024];
	char hap[1024];
	sprintf(bed, "%s/%s.bed", argv[5], argv[6]);
    sprintf(hap, "%s/%s.fasta", argv[5], argv[6]);

	FILE *bed_f = fopen(bed, "w");
	if (!bed_f) {
		fprintf(stderr, "Could not open bed file!\n");
		exit(EXIT_FAILURE);
	}

    FILE *hap_f = fopen(hap, "w");
	if (!hap_f) {
		fprintf(stderr, "Could not open haplotype file!\n");
		exit(EXIT_FAILURE);
	}
	struct sv_t *svs = NULL;
	int n_sv = 0;

	fprintf(stderr, "Loading genome ... \n");
	struct genome_t genome = read_genome(genome_fp);

	/* debug load genome */
	// fprintf(stderr, "%d\n", genome.sz);
	// for (int i = 0; i < genome.sz; ++i)
	// 	fprintf(stderr, "%s\t%d\n", genome.ref_name[i], genome.chr_sz[i]);

	fprintf(stderr, "Parsing variant ... \n");
	read_retro(retro_fp, &genome);
	read_sv(sv_fp, &svs, &n_sv, &genome);
	read_snp(snp_fp, &svs, &n_sv, &genome);

	fprintf(stderr, "Creating haplotype ... \n");
	sim_sv(svs, n_sv, &genome, hap_f, bed_f);

	return 0;
}
