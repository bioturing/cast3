#include "genome.h"
#include <assert.h>

inline std::string get_chr_name(std::string &s)
{
    return s.substr(1, (s.find(" ") == std::string::npos ? s.size() : s.find(" ")) - 1);
}

genome_t::genome_t(const char *fa_path, const char *bed_path, const char *fai_path)
{
    fprintf(stderr, "Loading genome ... ");
    std::ifstream fi(fa_path);
    if (!fi.is_open()) {
        perror("\nCould not open genome file");
        exit(EXIT_FAILURE);
    }
    std::string s, seq, cur_name;
    int64_t pos = 0;
    chr.clear();
    chr_pos.clear();
    ref_name.clear();

    while (getline(fi, s)) {
        if (s[0] == '>') {
            if (cur_name.empty()) {
                cur_name = get_chr_name(s);
                seq.clear();
                continue;
            }
            chr.push_back(seq);
            pos += seq.size();
            chr_pos.push_back(pos);
            ref_name.push_back(cur_name);
            cur_name = get_chr_name(s);
            seq.clear();
        } else {
            seq += s;
        }
    }

    ref_name.push_back(cur_name);
    chr.push_back(seq);
    pos += seq.size();
    chr_pos.push_back(pos);
    segment.resize(chr_pos.size());
    fprintf(stderr, "done\n");
    fi.close();

    fi.open(fai_path);
    if (!fi.is_open()) {
        perror("\nCould not open index reference file");
        exit(EXIT_FAILURE);
    }

    int64_t temp;
    for (int i = 0; i < (int)chr_pos.size(); ++i) {
        fi >> cur_name >> pos >> temp >> temp >> temp;
        assert(cur_name == ref_name[i]);
        segment[i].push_back(segment_t(1, pos + 1, cur_name, NOR));
    }

    load_bed(bed_path);
}

std::pair<int, int> genome_t::check_pos(int64_t pos, int len)
{
    auto up = std::upper_bound(chr_pos.begin(), chr_pos.end(), pos);
    if (pos + len >= *up)
        return std::make_pair(-1, -1);
    int tid = up - chr_pos.begin();
    return std::make_pair(tid, tid == 0 ? pos : pos - *(--up));
}

int64_t genome_t::get_chr_pos_back()
{
    return chr_pos.back();
}

std::string genome_t::get_ref_name(int id)
{
    assert(id < (int)ref_name.size());
    return ref_name[id];
}

std::string genome_t::get_sub_chr(int id, int pos, int len)
{
    assert(id < (int)ref_name.size());
    assert(pos + len <= (int)chr[id].size());
    return chr[id].substr(pos, len);
}

std::string get_type(int type)
{
    switch (type) {
    case NOR: return "NOR"; break;
    case TRA: return "TRA"; break;  
    case INV: return "INV"; break;
    case DEL: return "DEL"; break;
    case INS: return "INS"; break;
    case DUP: return "DUP"; break;
    }
    assert(false);
}

void process_inv(const std::string &ref, int lpos, int rpos, std::vector<segment_t> &chr)
{
    int found = 0;
    int old_start, old_end;
    std::vector<segment_t>::iterator it;

    for (it = chr.begin(); it != chr.end(); ++it) {
        if (it->type != NOR)
            continue;
        if (it->start <= lpos && it->end >= rpos) {
            old_start = it->start;
            old_end = it->end;
            it = chr.erase(it);
            it = chr.insert(it, segment_t(rpos, old_end, ref, NOR));
            it = chr.insert(it, segment_t(rpos - 1, lpos - 1, ref, INV));
            chr.insert(it, segment_t(old_start, lpos, ref, NOR));
            found = 1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr, "%s %d %d\n", ref.c_str(), lpos, rpos);
        fprintf(stderr, "Segment not found (INV)!\n");
        exit(EXIT_FAILURE);
    }
}

void process_del(const std::string &ref, int lpos, int rpos, std::vector<segment_t> &chr)
{
    int found = 0;
    int old_start, old_end;
    std::vector<segment_t>::iterator it;

    for (it = chr.begin(); it != chr.end(); ++it) {
        if (it->type != NOR)
            continue;
        if (it->start <= lpos && it->end >= rpos) {
            old_start = it->start;
            old_end = it->end;
            it = chr.erase(it);
            it = chr.insert(it, segment_t(rpos, old_end, ref, NOR));
            chr.insert(it, segment_t(old_start, lpos, ref, NOR));
            found = 1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr, "%s %d %d\n", ref.c_str(), lpos, rpos);
        fprintf(stderr, "Segment not found (DEL)!\n");
        exit(EXIT_FAILURE);
    }
}

void process_tra(const std::string &ref, int lpos, int rpos, const std::string &n_ref,
                 int n_lpos, int n_rpos, std::vector<segment_t> &chr)
{
    int found = 0;
    int old_start, old_end;
    std::vector<segment_t>::iterator it;

    for (it = chr.begin(); it != chr.end(); ++it) {
        if (it->type != NOR)
            continue;
        if (it->start <= lpos && it->end >= rpos) {
            old_start = it->start;
            old_end = it->end;
            it = chr.erase(it);
            it = chr.insert(it, segment_t(rpos, old_end, ref, NOR));
            it = chr.insert(it, segment_t(n_lpos, n_rpos, n_ref, TRA));
            chr.insert(it, segment_t(old_start, lpos, ref, NOR));
            found = 1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr, "%s %d %d\n", ref.c_str(), lpos, rpos);
        fprintf(stderr, "Segment not found (TRA)!\n");
        exit(EXIT_FAILURE);
    }
}

void process_ins(const std::string &ref, int lpos, int rpos, std::vector<segment_t> &chr)
{
    int found = 0;
    int old_start, old_end;
    std::vector<segment_t>::iterator it;

    for (it = chr.begin(); it != chr.end(); ++it) {
        if (it->type != NOR)
            continue;
        if (it->start <= lpos && it->end >= rpos) {
            old_start = it->start;
            old_end = it->end;
            it = chr.erase(it);
            it = chr.insert(it, segment_t(lpos, old_end, ref, NOR));
            it = chr.insert(it, segment_t(0, rpos - lpos, ref, INS));
            chr.insert(it, segment_t(old_start, lpos, ref, NOR));
            found = 1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr, "%s %d %d\n", ref.c_str(), lpos, rpos);
        fprintf(stderr, "Segment not found (INS)!\n");
        exit(EXIT_FAILURE);
    }
}

void process_dup(const std::string &ref, int lpos, int rpos, std::vector<segment_t> &chr)
{
    int found = 0;
    int old_start, old_end;
    std::vector<segment_t>::iterator it;

    for (it = chr.begin(); it != chr.end(); ++it) {
        if (it->type != NOR)
            continue;
        if (it->start <= lpos && it->end >= rpos) {
            old_start = it->start;
            old_end = it->end;
            it = chr.erase(it);
            it = chr.insert(it, segment_t(rpos, old_end, ref, NOR));
            it = chr.insert(it, segment_t(lpos, rpos, ref, DUP));
            chr.insert(it, segment_t(old_start, rpos, ref, NOR));
            found = 1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr, "%s %d %d\n", ref.c_str(), lpos, rpos);
        fprintf(stderr, "Segment not found (DUP)!\n");
        exit(EXIT_FAILURE);
    }
}

int genome_t::get_tid(const std::string &s)
{
    for (int i = 0; i < (int)ref_name.size(); ++i)
        if (ref_name[i] == s)
            return i;
    return -1;
}

void genome_t::load_bed(const char *file_path)
{
    fprintf(stderr, "Load bed ... ");
    std::ifstream fi;
    fi.open(file_path);
    if (!fi.is_open()) {
        perror("\nCould not open bed file");
        exit(EXIT_FAILURE);
    }

    std::string lref, rref, n_lref, n_rref, type, n_type;
    int lpos, rpos, ltid, rtid, n_lpos, n_rpos, i;

    while (fi >> lref >> lpos >> rref >> rpos >> type) {
        --lpos, --rpos;
        if (type == "SNP")
            continue;
        if (type == "INV") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_inv(lref, lpos, rpos, segment[ltid]);
            continue;
        }
        if (type == "DEL") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_del(lref, lpos, rpos, segment[ltid]);
            continue;
        }
        if (type == "TRA") {
            fi >> n_lref >> n_lpos >> n_rref >> n_rpos >> n_type;
            --n_lpos, --n_rpos;
            assert(n_type == type);
            assert(lref == n_lref && rref == n_rref);
            assert(lpos < n_lpos && rpos < n_rpos);
            ltid = get_tid(lref);
            rtid = get_tid(rref);
            if (ltid == rtid)
                assert(lpos > n_rpos || n_lpos < rpos);

            process_tra(lref, lpos, n_lpos, rref, rpos, n_rpos, segment[ltid]);
            process_tra(rref, rpos, n_rpos, lref, lpos, n_lpos, segment[rtid]);
        }
        if (type == "DUP") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_dup(lref, lpos, rpos, segment[ltid]);
            continue;
        }
        if (type == "INS") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_ins(lref, lpos, rpos, segment[ltid]);
            continue;
        }
    }

    for (i = 0; i < (int)segment.size(); ++i) {
        int pos = 0;
        for (auto &item : segment[i]) {
            item.pos = pos;
            pos += item.end - item.start;
        }
    }

    /* debug purpose */
    // for (int i = 0; i < (int)genome.size(); ++i) {
    //  cout << i << endl;
    //  for (auto &x : genome[i]) {
    //      cout << "\t" << x.start << ", " << x.end << ", " << x.ref
    //           << ", " << x.pos << ", " << get_type(x.type) << "\n";
    //  }
    // }
    // cout << "----------------------------------\n";

    fi.close();
    fprintf(stderr, "done\n");
}

std::string genome_t::convert(int tid, int pos1, int len1, int pos2, int len2)
{
    int lo, hi, mid, is_first;
    std::string new_pos1, new_pos2;
    auto &chr = segment[tid];

    /* first pos */
    lo = 0, hi = (int)chr.size() - 1;
    while (lo < hi) {
        mid = ((lo + hi) >> 1) + 1;
        if (chr[mid].pos <= pos1)
            lo = mid;
        else
            hi = mid - 1;
    }

    is_first = 1;
    for ( ; ; ++lo) {
        int diff = pos1 - chr[lo].pos;
        int len = chr[lo].end - chr[lo].start;
        if (is_first)
            is_first = 0;
        else
            new_pos1 += "_";
        if (len - diff >= len1) {
            new_pos1 += chr[lo].ref + ":" + std::to_string(chr[lo].start + diff) +
                        ":" + std::to_string(len1) + ":" + get_type(chr[lo].type);
            break;
        }
        new_pos1 += chr[lo].ref + ":" + std::to_string(chr[lo].start + diff) + ":" +
                    std::to_string(len - diff) + ":" + get_type(chr[lo].type);
        len1 -= len - diff;
        assert(lo + 1 < (int)chr.size());
        pos1 = chr[lo + 1].pos;
    }

    /* second pos */
    lo = 0, hi = (int)chr.size() - 1;
    while (lo < hi) {
        mid = ((lo + hi) >> 1) + 1;
        if (chr[mid].pos <= pos2)
            lo = mid;
        else
            hi = mid - 1;
    }

    is_first = 1;
    for ( ; ; ++lo) {
        int diff = pos2 - chr[lo].pos;
        int len = chr[lo].end - chr[lo].start;
        if (is_first)
            is_first = 0;
        else
            new_pos2 += "_";
        if (len - diff >= len2) {
            new_pos2 += chr[lo].ref + ":" + std::to_string(chr[lo].start + diff) +
                        ":" + std::to_string(len2) + ":" + get_type(chr[lo].type);
            break;
        }
        new_pos2 += chr[lo].ref + ":" + std::to_string(chr[lo].start + diff) + ":" +
                    std::to_string(len - diff) + ":" + get_type(chr[lo].type);
        len2 -= len - diff;
        assert(lo + 1 < (int)chr.size());
        pos2 = chr[lo + 1].pos;
    }

    /* result after convert */
    return "@" + new_pos1 + "|" + new_pos2;
}