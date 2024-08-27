#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "mio.hpp"

#define MAX_READ_SIZE 2500000
#ifndef PRINT_CHROSOMAL
#define PRINT_CHROSOMAL 0
#endif
using std::cout;
using std::ostream;
using std::set;
using std::string;
using std::vector;

int process_fastq(int argc, char **argv);
int process_gaf(int argc, char **argv);

int main(int argc, char **argv) { return process_gaf(argc, argv); }
struct gview {
    unsigned long s;
    unsigned long e;
};
struct mmap_view {
    const mio::mmap_source *mmap;
    unsigned long s;
    unsigned long e;
    unsigned long cap;

    mmap_view(const mio::mmap_source *mmap)
        : mmap{mmap}, s{0}, e{0}, cap{mmap->size()} {}
    mmap_view(const mio::mmap_source *mmap, unsigned long s, unsigned long e,
              unsigned long cap)
        : mmap{mmap}, s{s}, e{e}, cap{cap} {}
    mmap_view(const mio::mmap_source *mmap, gview g)
        : mmap{mmap}, s{g.s}, e{g.e} {}
    friend ostream &operator<<(ostream &ost, mmap_view &view) {
        for (unsigned long i = view.s; i < view.e; ++i) {
            ost << (char)(*view.mmap)[i];
        }
        return ost;
    }
    operator string() { return string{mmap->begin() + s, mmap->begin() + e}; }
    operator gview() { return {s, e}; }
    bool operator==(const kstring_t &other) const {
        if (other.l != e - s) {
            return false;
        }
        for (unsigned i = s; i < e; ++i) {
            unsigned j = i - s;
            if ((*mmap)[i] != other.s[j]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const string &other) const {
        if (other.size() != e - s) {
            return false;
        }
        for (unsigned i = s; i < e; ++i) {
            unsigned j = i - s;
            if ((*mmap)[i] != other[j]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const mmap_view &other) const {
        if (other.e - other.s != e - s) {
            return false;
        }
        for (unsigned i = s; i < e; ++i) {
            unsigned j = i - s;
            if ((*mmap)[i] != (*other.mmap)[j + other.s]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const gview &other) const {
        if (other.e - other.s != e - s) {
            return false;
        }
        for (unsigned i = s; i < e; ++i) {
            unsigned j = i - s;
            if ((*mmap)[i] != (*mmap)[j + other.s]) {
                return false;
            }
        }
        return true;
    }
    mmap_view focus() { return mmap_view{mmap, s, s, e}; }

    void skip_next(string delimiters) {
        while (s < cap) {
            for (char delimiter : delimiters) {
                if ((*mmap)[s] == delimiter) {
                    s++;
                    e = s;
                    return;
                }
            }
            ++s;
        }
        s++;
        e = s;
    }
    void skip_next(char delimiter) {
        while (s < cap) {
            if ((*mmap)[s] == delimiter) {
                break;
            }
            ++s;
        }
        s++;
        e = s;
    }

    char extend_until(string delimiters) {
        while (e < cap) {
            for (char delimiter : delimiters) {
                if ((*mmap)[e] == delimiter) {
                    return delimiter;
                }
            }
            ++e;
        }
        return 0;
    }
    void extend_until(char delimiter) {
        while (e < cap) {
            if ((*mmap)[e] == delimiter) {
                break;
            }
            ++e;
        }
    }
    void skip_next_n(char delimiter, int N) {
        for (int i = 0; i < N; ++i) {
            this->skip_next(delimiter);
        }
    }
};


set<string> parse_plasmid_tsv(const string &path) {
    set<string> plasmid_contigs;
    std::error_code error;
    mio::mmap_source tsv_mmap = mio::make_mmap_source(path, error);
    if (error) {
        const auto &errmsg = error.message();
        std::fprintf(stderr, "error mapping file: %s, exiting...\n",
                     errmsg.c_str());
        exit(-1);
    }

    mmap_view tsv_view{&tsv_mmap};
    tsv_view.skip_next('\n');  // Header
    while (tsv_view.e < tsv_view.cap) {
        tsv_view.extend_until("\t\n");

        plasmid_contigs.insert(tsv_view);
        tsv_view.skip_next('\n');
    }
    return plasmid_contigs;
}

struct GafEntry {
    bool isterm;
    gview id;
    vector<gview> contigs;
    GafEntry(mmap_view &gaf_view) : isterm{false}, id{(unsigned)-1,(unsigned)-1} {
        bool first = true;
        if (gaf_view.e >= gaf_view.cap) {
            isterm = true;
            return;
        }
        do {
            gaf_view.extend_until('\t');
            if (gaf_view != id) {
                if (!first) {
                    gaf_view.e = gaf_view.s;  // rollback
                    break;
                }
                id = gaf_view;
            }
            first = false;
            gaf_view.skip_next_n('\t', 5);
            gaf_view.extend_until('\t');
            mmap_view alig_view = gaf_view.focus();
            alig_view.skip_next("<>");

            while (alig_view.e < alig_view.cap) {
                alig_view.extend_until("<>");
                contigs.push_back(alig_view);
                alig_view.skip_next("<>");
            }

            gaf_view.skip_next('\n');
        } while (gaf_view.e < gaf_view.cap);
    }
};
int process_gaf(int argc, char **argv) {
    std::error_code error;

    mio::mmap_source gaf_mmap = mio::make_mmap_source(argv[1], error);
    if (error) {
        const auto &errmsg = error.message();
        std::fprintf(stderr, "error mapping file: %s, exiting...\n",
                     errmsg.c_str());
        return error.value();
    }

    gzFile fastq_fp = gzopen(argv[2], "r");
//    gzFile plasmid_out = gzopen(argv[4], "w");

    FILE *plasmid_out = popen((string{"gzip - > "} + string{argv[4]}).c_str(), "w");
    FILE *unknown_out = popen((string{"gzip - > "} + string{argv[6]}).c_str(), "w");

#if PRINT_CHROSOMAL
    FILE *chrmsm_out = popen((string{"gzip - > "} + string{argv[5]}).c_str(), "w");
#endif
    set<string> plasmid_contigs = parse_plasmid_tsv(argv[3]);
    kseq_t *seq = kseq_init(fastq_fp);

    mmap_view gaf_view{&gaf_mmap};

    string prev_id;

    char *buffer = new char[MAX_READ_SIZE];

    int l;
    int buffer_length = 0;


    //    l = kseq_read(seq);

    vector<string> current_contigs;

    for (GafEntry g{gaf_view}; !g.isterm; g = GafEntry{gaf_view}) {
        l = kseq_read(seq);
        mmap_view gid {&gaf_mmap, g.id};
        while (l >= 0) {

            if (gid != seq->name) {
                buffer_length = fprintf(unknown_out, "@%s %s\n%s\n+\n%s\n", seq->name.s,
                                    seq->comment.s, seq->seq.s, seq->qual.s);


                l = kseq_read(seq);
            } else {
                break;
            }
        }
        bool from_chromosome = true;
        for (const gview v : g.contigs) {
            mmap_view vv{&gaf_mmap, v};
            if (plasmid_contigs.contains((string)vv)) {
                buffer_length = fprintf(plasmid_out, "@%s %s\n%s\n+\n%s\n", seq->name.s,
                                    seq->comment.s, seq->seq.s, seq->qual.s);

                from_chromosome = false;
                break;
            }
        }
#if PRINT_CHROSOMAL
        if(from_chromosome){
            buffer_length = fprintf(chrmsm_out, "@%s %s\n%s\n+\n%s\n", seq->name.s,
                                seq->comment.s, seq->seq.s, seq->qual.s);


        }
#endif
    }

    kseq_destroy(seq);
    gzclose(fastq_fp);


    pclose(plasmid_out);
    pclose(unknown_out);
#if PRINT_CHROSOMAL
    pclose(chrmsm_out);
#endif
    delete[] buffer;
    return 0;
}

int process_fastq(int argc, char **argv) {
    //    FILE *gaf_fp  = fopen(argv[1], "r");
    gzFile fastq_fp = gzopen(argv[1], "r");
    gzFile plasmid_out = gzopen(argv[2], "w");
    //    gzFile chr_out = gzopen(argv[4], "w");

    set<string> plasmid_contigs;

    kseq_t *seq = kseq_init(fastq_fp);
    int l;
    l = kseq_read(seq);
    int cnt = 0;
    char buffer[2500000];

    while (l >= 0) {
        int ll = sprintf(buffer, "@%s %s\n%s\n+\n%s\n", seq->name.s,
                         seq->comment.s, seq->seq.s, seq->qual.s);
        gzwrite(plasmid_out, buffer, ll);
        //        int wr=gzprintf(plasmid_out, "@%s %s\n%s\n+\n%s\n",
        //        seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);

        l = kseq_read(seq);
        ++cnt;
    }

    printf("return value: %d\n%d reads\n", l, cnt);
    kseq_destroy(seq);
    gzclose(fastq_fp);
    gzclose(plasmid_out);
    return 0;
}
