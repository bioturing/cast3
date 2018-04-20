#ifndef _GENERATE_H_
#define _GENERATE_H_

#include <random>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <algorithm>

#include "genome.h"

#define BARCODE_SZ                  16

#define MEAN_MLC_COVERAGE           20

#define GAMMA_DIS_A                 2.178355
#define GAMMA_DIS_B                 1 / 2.541399e-5

#define NORMAL_MEAN                 351
#define NORMAL_STDDEV               50

#define READ1_ERR                   0.005
#define READ2_ERR                   0.005

#define MIN_MLC_LEN                 1000
#define MAX_MLC_LEN                 1000000

struct molecule_t {
    int n_read_pair;
    int len;
    int start_pos;
    int tid;
    int hap;

    molecule_t(int _n_read_pair, int _len)
    {
        n_read_pair = _n_read_pair;
        len = _len;
        start_pos = -1;
        hap = -1;
        tid = -1;
    }
};

class generate_t {
private:
    std::vector<std::string> barcode_lst;
    genome_t hap1, hap2;
    std::ofstream fi_mlc;

    molecule_t generate_molecule(int read_len);
    std::vector<molecule_t> generate_barcode(int read_len, int &total_read, int mean_mlc_per_bx);

public:
    generate_t(const char *bx_path, const char *hap1_path, const char *bed1_path,
               const char *hap2_path, const char *bed2_path, const char *fai_path);
    void generate_read(int read_len, int total_read, int mean_mlc_per_bx);
};

#endif /* _GENERATE_H_ */
