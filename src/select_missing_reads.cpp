#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include <iostream>
#include <set>
#include <unordered_map>

#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "mio.hpp"

#define DEBUG 1

#ifdef DEBUG
#define LOGf(...) fprintf(stderr, __VA_ARGS__)
#else
#define LOGf(...) 
#endif
#define MAX_READ_SIZE 2500000
#ifndef PRINT_CHROSOMAL
#define PRINT_CHROSOMAL 0
#endif
using std::ostream;
using std::set;
using std::unordered_map;
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

    char operator [] (size_t idx){
        return (*mmap)[idx + s];
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

    char at_end(){
        if ( e == cap){
            return 0;
        }
        else{
            return (*mmap)[e];
        }
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

enum class ContigType {
    Plasmid,
    Chromosome,
    Unknown
};

unordered_map<string, ContigType> parse_plasmid_tsv(const string &path){
    unordered_map<string, ContigType> plasmid_contigs;
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
        char nextc = tsv_view.at_end();
        ContigType ctgt = ContigType::Plasmid;
        string ctg = tsv_view;
        if (nextc == '\t'){
            tsv_view.skip_next('\t');
            tsv_view.extend_until("\t\n");
            if (tsv_view == "plasmid"){
                ctgt = ContigType::Plasmid;
            }
            else if (tsv_view == "chromosome"){
                ctgt = ContigType::Chromosome;
            }
            else{
                ctgt = ContigType::Unknown;
            }
        }
        LOGf("%s\t.\t%s\n",ctg.c_str(), ((string) tsv_view).c_str());
        plasmid_contigs[ctg] = ctgt;
        tsv_view.skip_next('\n');
    }
    return plasmid_contigs;
}

set<string> parse_plasmid_tsv_old(const string &path) {
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

struct PafEntry {
    bool isterm;
    gview id;
    vector<gview> contigs;
    PafEntry(mmap_view &paf_view) : isterm{false}, id{(unsigned)-1,(unsigned)-1} {
        bool first = true;
        if (paf_view.e >= paf_view.cap) {
            isterm = true;
            return;
        }
        do {
            paf_view.extend_until('\t');
            if (paf_view != id) {
                if (!first) {
                    paf_view.e = paf_view.s;  // rollback
                    break;
                }
                id = paf_view;
            }
            first = false;
            paf_view.skip_next_n('\t', 5);
            paf_view.extend_until('\t');
            mmap_view alig_view = paf_view.focus();

            if (alig_view[0] == '<' || alig_view[0] == '>'){
                alig_view.skip_next("<>");
            }
            while (alig_view.e < alig_view.cap) {
                alig_view.extend_until("<>");
                contigs.push_back(alig_view);
                alig_view.skip_next("<>");
            }

            paf_view.skip_next('\n');
        } while (paf_view.e < paf_view.cap);
    }
};
int process_gaf(int argc, char **argv) {
    std::error_code error;

    mio::mmap_source paf_mmap = mio::make_mmap_source(argv[1], error);
    if (error) {
        const auto &errmsg = error.message();
        std::fprintf(stderr, "error mapping file: %s, exiting...\n",
                     errmsg.c_str());
        return error.value();
    }
    mio::mmap_source gaf_mmap = mio::make_mmap_source(argv[2], error);
    if (error) {
        const auto &errmsg = error.message();
        std::fprintf(stderr, "error mapping file: %s, exiting...\n",
                     errmsg.c_str());
        return error.value();
    }


    gzFile fastq_fp = gzopen(argv[3], "r");

    FILE *plasmid_out = popen((string{"gzip - > "} + string{argv[4]}).c_str(), "w");

    auto blacklist_contigs = parse_plasmid_tsv(argv[5]);
    

    kseq_t *seq = kseq_init(fastq_fp);


    mmap_view paf_view{&paf_mmap};
    mmap_view gaf_view{&gaf_mmap};

    string prev_id;

    char *buffer = new char[MAX_READ_SIZE];

    //    l = kseq_read(seq);

    vector<string> current_contigs;
    unordered_map<string, int> reads2use;

    for (PafEntry g{paf_view}; !g.isterm; g = PafEntry{paf_view}) {
        mmap_view gid {&paf_mmap, g.id};

        for (const gview v : g.contigs) {
            mmap_view vv{&paf_mmap, v};
            reads2use[(string)vv] = 1;
        }
    }

    for (PafEntry g{gaf_view}; !g.isterm; g = PafEntry{gaf_view}) {
        mmap_view gid {&gaf_mmap, g.id};

        for (const gview v : g.contigs) {
            mmap_view vv{&gaf_mmap, v};
            auto it = blacklist_contigs.find((string) vv);
            if(it != blacklist_contigs.end() && it->second == ContigType::Chromosome){
                auto it2 = reads2use.find((string)vv);
                if(it2!=reads2use.end()){
                    it2->second = 0;
                }
            }
        }
    }


    while (kseq_read(seq) >= 0){

        string name{seq->name.s};
        auto it = reads2use.find(name);
        if(it != reads2use.end() && it->second > 0){
            fprintf(plasmid_out, "@%s %s\n%s\n+\n%s\n", seq->name.s,
                                seq->comment.s, seq->seq.s, seq->qual.s);
            it->second--;
        }
    }
    kseq_destroy(seq);
    gzclose(fastq_fp);


    pclose(plasmid_out);


    delete[] buffer;
    return 0;
}

