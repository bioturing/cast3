#ifndef _GENOME_H_
#define _GENOME_H_

#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

#define NOR                         0
#define INV                         1
#define INS                         2
#define DUP                         3
#define DEL                         4

struct segment_t {
    /* start - end is 1-base, pos is 0-base */
    int start;
    int end;
    std::string ref;
    int type;
    int pos;

    segment_t(int _start, int _end, const std::string &_ref, int _type)
    {
        start = _start;
        end = _end;
        ref = _ref;
        type = _type;
        pos = -1;
    }
};

class genome_t {
private:
    std::vector<std::string> chr;
    std::vector<std::string> ref_name;
    std::vector<int64_t> chr_pos;
    std::vector<std::vector<segment_t> > segment;
    int id;

    void load_bed(const char *file_path);

public:
    genome_t(const char *fa_path, const char *bed_path, const char *fai_path, int _id);
    std::pair<int, int> check_pos(int64_t pos, int len);
    int64_t get_chr_pos_back();
    int get_tid(const std::string &s);
    std::string get_ref_name(int id);
    int get_chr_size(int id);
    std::string get_sub_chr(int id, int pos, int len);
    std::string convert(int tid, int pos1, int len1, int pos2, int len2);
};

#endif /* _GENOME_H_ */