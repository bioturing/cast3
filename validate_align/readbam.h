#ifndef _READBAM_H_
#define _READBAM_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "htslib/sam.h"

int get_tid(char *s);

char *get_name(int tid);

void sam_write(FILE *fp, bam1_t *b);

char *convert_scigar(uint32_t *bcigar, int sz);

void init_bam(char *file_path);

int read_bam(bam1_t *b);

#endif /* _READBAM_H_ */