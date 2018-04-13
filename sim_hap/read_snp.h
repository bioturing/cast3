#ifndef _READ_SNP_H_
#define _READ_SNP_H_ 

#include <stdio.h>
#include <stdlib.h>

struct sv_t;

/* read SNP from bed file then process SNPs and store indels
 * as the structural variants
 */
void read_snp(FILE *fp, struct sv_t *svs, int *n, char **genome, FILE *bed_f);

#endif /* _READ_SNP_H_ */
