#ifndef READ_SV_H
#define READ_SV_H

#define MAX_SV			3000000
#define CNV_POS			3
#define MAX_DUP 		3
#define BUFFER_SIZE 		1024
#define NAME_BUFFER 		100
#define MAX_RETROS 		3000

#define GET_CNV_(alt) alt[CNV_POS] - '0'
#define GET_CNV(alt) GET_CNV_(alt) > MAX_DUP ? MAX_DUP : GET_CNV_(alt)

#define IS_DEL(sv) (!strcmp(sv->type, "CNV")) && sv->cnv == 0
#define IS_DUP(sv) (!strcmp(sv->type, "CNV")) && sv->cnv > 0
#define IS_INV(sv) (!strcmp(sv->type, "INV"))
#define IS_INS(sv) (!strcmp(sv->type, "INS"))
#define IS_SINS(sv) (!strcmp(sv->type, "SINS"))
#define IS_RETRO(sv) !strcmp(sv->type, "ALU") || !strcmp(sv->type, "SVA") || \
	!strcmp(sv->type, "LINE1")

#define CHROM_IDX(chr, i) for(i = 0; strcmp(chr, *(chroms + i)); ++i)

#include "read_snp.h"
#include "read_fasta.h"

const uint8_t nst_nt4_table[256];

/* structural variant struct
 * @att chr: chromosome idx in the chroms array
 * @att type: type of event CNV, INS, DEL = CNV with cnv=0, DUP = CNV with cnv>0
 * @att start: start position of the event
 * @att end: end position of the event
 * @att cnv: copy number of the event, only use for DEL and DUP
 * @att seq: seq
 */
struct sv_t{
	int chr;
	char type[20];
	int start;
	int end;
	int cnv;
	char *seq;
};

struct seg_t{
	int start;
	int end;
	char *seq;
	char strand;
};

/* Read the structural variants from the bed file fp, save the number of structural 
 * variants
 * @param fp: bed file
 * @param 
void read_sv(FILE *fp, struct sv_t *svs, int *n);

 * Simple normal random number generator, copied from genran.c 
 * @return random double in range [0,1]
 * */
double ran_normal();

/* Write the sequence of each segment chain into the out stream
 * @param seg: an array of segments
 * @param idx: index of segment chain,e.i index of the corresponding chromosome
 * @param n: length of seg
 * @param genome: the genome sequences
 * @param fp: out stream
 * @return: void, write sequence into the out stream in fasta format
 */
void seg2fasta(struct seg_t *seg, int idx, int n, char **genome, FILE *fp);

/* Sum up each segment size in the chain of segment
 * @param seg: an array of segments
 * @param n: length of seg
 * @return : total size of segments
 */
int segsize(struct seg_t *seg, int n);

/* Compare each structural variant by their chromosome first, then by their start
 * positions
 * @return: integer value determines whether a < b, a = b, or a > b
 */
int comp_sv(const void *a, const void *b);

/* Main function for generate haplotype with SNPs, retrotransposon, and structural
 * variants.
 * Usage: 1. Read the genome, func: read_genome
 * 		  2. Read SV, SNPs, and retro transposons, func: read_sv, read_snp, read_retro
 * 		  3. Simulate the haplotype, func: sim_sv
 * @param svs: array of structural variants
 * @param n: length of svs
 * @param seq: the genome
 * @param segs: segments for all chromosomes in the genome, each chrom has an
 * 				segment chain
 * @param n_segs: array length of each segment chain
 * @param stream: out stream of the haplotype sequence
 */
void sim_sv(struct sv_t *svs, int n, char **seq, struct seg_t **segs, int *n_segs,
				 FILE *stream, FILE *bed_f);
void read_retro(FILE *fp, char **seqs);

#endif /* READ_SV_H */
