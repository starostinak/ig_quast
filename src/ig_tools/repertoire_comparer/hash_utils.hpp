#pragma once

#include "../utils/string_tools.hpp"
#include <unordered_map>

struct kmer_match {
    size_t id;
    string const * seq;
    size_t pos;
    bool strand;
};

typedef unordered_map <size_t, vector <kmer_match> > KmerHashMap;

class HashBuilder {
public:
    typedef size_t Hash;
    typedef size_t Pos;
    HashBuilder(size_t hash_size, size_t indent) : hash_size_(hash_size), indent_(indent) { }

    vector <pair <Pos, Hash> > get_hashes(string const & seq) const {
        vector <pair <Pos, Hash> > res_hashes;
        vector <size_t> seq_hashes;
        StringToKHashes(seq, hash_size_, seq_hashes);
        for (size_t i = 0; i != seq_hashes.size(); ++i) {
            bool use_hash = true;
            if (indent_ && i >= indent_ - hash_size_ && i < seq.size() - indent_) {
                use_hash = false;
            }
            if (use_hash) {
                res_hashes.push_back(make_pair(i, seq_hashes[i]));
            }
        }
        return res_hashes;
    }

private:
    size_t hash_size_;
    size_t indent_;
};

void build_kmers_hash_map(map <size_t, string> const & seqs, HashBuilder const & hash_builder,
                          KmerHashMap & kmers_hash_map, bool strand) {
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        vector <pair <size_t, size_t> > pos_hashes = hash_builder.get_hashes(it->second);
        for (auto pos_hash_it = pos_hashes.begin(); pos_hash_it != pos_hashes.end(); ++pos_hash_it) {
            kmers_hash_map[pos_hash_it->second].push_back({ it->first, &it->second,
                                                            pos_hash_it->first, strand });
        }
    }
}
