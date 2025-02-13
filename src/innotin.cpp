#include <iostream>
#include <unordered_set>
#include <zlib.h>
#include "kseq.h"

// Initialize kseq for reading gzipped or plain FASTQ files
KSEQ_INIT(gzFile, gzread)

// Function to read FASTQ files and extract read IDs
std::unordered_set<std::string> load_read_ids(const std::string& filename) {
    std::unordered_set<std::string> read_ids;
    gzFile fp = gzopen(filename.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Unable to open file " << filename << "\n";
        exit(1);
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        read_ids.insert(seq->name.s);  // Store read ID
    }

    kseq_destroy(seq);
    gzclose(fp);
    return read_ids;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " main.fastq subset1.fastq [subset2.fastq ...]\n";
        return 1;
    }

    std::string main_fastq = argv[1];

    // Load read IDs from all subset FASTQs into a set
    std::unordered_set<std::string> subset_reads;
    for (int i = 2; i < argc; ++i) {
        std::unordered_set<std::string> reads = load_read_ids(argv[i]);
        subset_reads.insert(reads.begin(), reads.end());
    }

    // Process the main FASTQ file and print only the unique reads
    gzFile fp = gzopen(main_fastq.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: Unable to open main FASTQ file " << main_fastq << "\n";
        return 1;
    }

    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        if (subset_reads.find(seq->name.s) == subset_reads.end()) {
            std::cout << ">" << seq->name.s << " " << seq->comment.s << "\n"
                      << seq->seq.s << "\n";
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    return 0;
}

