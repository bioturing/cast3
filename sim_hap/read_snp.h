#ifndef _READ_SNP_H_
#define _READ_SNP_H_ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "read_sv.h"
#include "read_fasta.h"

/* read SNP from bed file then process SNPs and store indels
 * as the structural variants
 */
void read_snp(FILE *fp, struct sv_t **svs, int *n_sv, struct genome_t *genome);

#endif /* _READ_SNP_H_ */
