#ifndef _READ_SV_H_
#define _READ_SV_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "read_fasta.h"

#define IS_DEL(sv)		((strcmp(sv->type, "CNV") == 0) && sv->cnv == 0)
#define IS_DUP(sv)		((strcmp(sv->type, "CNV") == 0) && sv->cnv > 0)
#define IS_INV(sv)		(strcmp(sv->type, "INV") == 0)
#define IS_INS(sv)		(strcmp(sv->type, "INS") == 0)
#define IS_SINS(sv)		(strcmp(sv->type, "SINS") == 0)
#define IS_RETRO(sv)		(strcmp(sv->type, "ALU") == 0 ||	\
				 strcmp(sv->type, "SVA") == 0 ||	\
				 strcmp(sv->type, "LINE1") == 0)

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
 */
void read_sv(FILE *fp, struct sv_t **svs, int *n_sv, struct genome_t *genome);

/* Simple normal random number generator, copied from genran.c 
 * @return random double in range [0,1]
 */
double ran_normal();

/* Write the sequence of each segment chain into the out stream
 * @param seg: an array of segments
 * @param n_seg: size of seg
 * @param idx: index of segment chain, e.i index of the corresponding chromosome
 * @param genome: the genome data
 * @param fp: out stream
 * @return: void, write sequence into the out stream in fasta format
 */
void seg2fasta(struct seg_t *seg, int n_seg, int idx, struct genome_t *genome, FILE *fp);

/* Sum up each segment size in the chain of segment
 * @param seg: an array of segments
 * @param n: size of seg
 * @return : total size of segments
 */
int segsize(struct seg_t *seg, int n);

/* Main function for generate haplotype with SNPs, retrotransposon, and structural
 * variants.
 * Usage: 1. Read the genome, func: read_genome
 * 	  2. Read SV, SNPs, and retro transposons, func: read_sv, read_snp, read_retro
 * 	  3. Simulate the haplotype, func: sim_sv
 * @param svs: array of structural variants
 * @param n_sv: size of svs
 * @param genome: the genome data
 * @param hap_f: out stream of the haplotype fasta
 * @param bed_f: out stream of the bed data
 */
void sim_sv(struct sv_t *svs, int n_sv, struct genome_t *genome, FILE *hap_f, FILE *bed_f);

void read_retro(FILE *fp, struct genome_t *genome);

#endif /* READ_SV_H */
