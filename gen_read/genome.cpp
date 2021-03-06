#include "genome.h"
#include <assert.h>

std::string get_chr_name(std::string &s)
{
    int pos = s.find(" ");
    if (pos == (int)std::string::npos)
        pos = s.find("\t");
    if (pos == (int)std::string::npos)
        pos = s.size();
    return s.substr(1, (s.find(" ") == std::string::npos ? s.size() : s.find(" ")) - 1);
}

genome_t::genome_t(const char *fa_path, const char *bed_path, const char *fai_path, int _id) : id(_id)
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
    assert(up != chr_pos.end());
    if (pos + len > *up)
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

int genome_t::get_chr_size(int id)
{
    assert(id < (int)chr.size());
    return chr[id].size();
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
            if (lpos < old_end)
                it = chr.insert(it, segment_t(lpos, old_end, ref, NOR));
            it = chr.insert(it, segment_t(lpos, rpos, ref, DUP));
            if (old_start < lpos)
                chr.insert(it, segment_t(old_start, lpos, ref, NOR));
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

    std::string lref, rref, type;
    int lpos, rpos, ltid, i;

    while (fi >> lref >> lpos >> rref >> rpos >> type) {
        ++lpos, ++rpos;
        if (type == "SNP") {
            continue;
        } else if (type == "INV") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_inv(lref, lpos, rpos, segment[ltid]);
            continue;
        } else if (type == "DEL") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_del(lref, lpos, rpos, segment[ltid]);
            continue;
        } else if (type == "DUP") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_dup(lref, lpos, rpos, segment[ltid]);
            continue;
        } else if (type == "INS") {
            assert(lref == rref);
            assert(lpos < rpos);
            ltid = get_tid(lref);
            process_ins(lref, lpos, rpos, segment[ltid]);
            continue;
        } else {
            assert(false);
        }
    }

    for (i = 0; i < (int)segment.size(); ++i) {
        int pos = 0;
        for (auto &item : segment[i]) {
            item.pos = pos;
            pos += abs(item.end - item.start);
        }
        // std::cerr << get_ref_name(i) << " " << pos << " " << (int)chr[i].size() << std::endl;
        assert(pos == (int)chr[i].size());
    }

    /* debug purpose */
    // std::ofstream fo("debug_segment");
    // for (int i = 0; i < (int)segment.size(); ++i) {
    //     fo << i << std::endl;
    //     for (auto &x : segment[i]) {
    //         fo << "\t" << x.start << ", " << x.end << ", " << x.ref
    //            << ", " << x.pos << ", " << get_type(x.type) << "\n";
    //     }
    // }
    // fo << "----------------------------------\n";
    // fo.close();

    fi.close();
    fprintf(stderr, "done\n");
}

std::string genome_t::convert(int tid, int pos1, int len1, int pos2, int len2)
{
    int lo, hi, mid, is_first;
    std::string new_pos1, new_pos2;
    auto &chr = segment[tid];

    int old_pos2 = pos2, old_pos1 = pos1;
    int old_len2 = len2, old_len1 = len1;

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
        if (lo + 1 == (int)chr.size())
            std::cerr << old_pos1 << " " << old_len1 << " " << get_ref_name(tid) << " "
                      << id << " " << std::endl;
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
        if (lo + 1 == (int)chr.size())
            std::cerr << old_pos2 << " " << old_len2 << " " << get_ref_name(tid) << " "
                      << id << " " << std::endl;
        assert(lo + 1 < (int)chr.size());
        pos2 = chr[lo + 1].pos;
    }

    /* result after convert */
    return new_pos1 + "|" + new_pos2;
}
