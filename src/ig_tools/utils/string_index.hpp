#pragma once

#include "include_me.hpp"
#include "string_tools.hpp"
#include <unordered_map>
#include <unordered_set>

class IndexSelector {
public:
    virtual bool NeedIncludeToIndex(const string & seq, size_t offset, size_t length) const = 0;
};

struct ShiftedString {
    bool operator == (const ShiftedString & other) const
    {
       return seq_number == other.seq_number && relative_shift == other.relative_shift; 
    }

    struct Hash {
        size_t operator () (const ShiftedString & other) const
        {
            return other.seq_number * 997 + other.relative_shift;
        }
    };

    size_t seq_number;
    int relative_shift;
};

class IndexedStringSet {
// TODO need new class name!
// strings are shifted and have errors.
public:
	IndexedStringSet(const vector <string> & sequences, // TODO is need to add some meta information?
                     const IndexSelector & index_selector,
                     size_t k_mer_length,
                     size_t min_overlap,
			         size_t max_k_mer_occurrences);

    void FindCandidates(const string & seq, vector <ShiftedString> & candidates) const;

private:
    static size_t GetSequencesOverlap(const string & seq1, const string & seq2, int relative_shift);

    struct KMerMatch {
        size_t seq_number;
        size_t offset;
    };

    const vector <string> & sequences;
    size_t k_mer_length;
    size_t min_overlap;
    size_t max_k_mer_occurrences;
    unordered_map < size_t, vector <KMerMatch> > matches_by_hash;
};

IndexedStringSet::IndexedStringSet(const vector <string> & sequences,
                                   const IndexSelector & index_selector,
                                   size_t k_mer_length,
                                   size_t min_overlap,
			                       size_t max_k_mer_occurrences)
    : sequences(sequences), k_mer_length(k_mer_length), min_overlap(min_overlap), max_k_mer_occurrences(max_k_mer_occurrences)
{
    for (size_t seq_number = 0; seq_number < sequences.size(); ++seq_number)
    {
        const string & seq = sequences[seq_number];
        vector <size_t> seq_hashes;
        StringToKHashes(seq, k_mer_length, seq_hashes);
        for (size_t offset = 0; offset < seq_hashes.size(); ++offset)
        {
            if (index_selector.NeedIncludeToIndex(seq, offset, k_mer_length))
            {
                matches_by_hash[seq_hashes[offset]].push_back({seq_number, offset});
            // TODO OPTIMIZATION: size of vector can exceed max_k_mer_occurences
            }
        }
    }
}

size_t IndexedStringSet::GetSequencesOverlap(const string & seq1, const string & seq2, int relative_shift)
{
    if (relative_shift >= 0)
    {
        return min(seq1.size() - relative_shift, seq2.size());
    }
    else
    {
        return min(seq1.size(), seq2.size() + relative_shift);
    }
}

void IndexedStringSet::FindCandidates(const string & seq, vector <ShiftedString> & candidates) const
{
    vector <size_t> seq_hashes;
    StringToKHashes(seq, k_mer_length, seq_hashes);
    size_t seq_length = seq.length();
    unordered_set < ShiftedString, ShiftedString::Hash > candidates_set;
    for (size_t i = 0; i < seq_length; ++i)
    {
        // TODO Skip some tricky Katya's check here
        auto matches_iterator = matches_by_hash.find(seq_hashes[i]);
        if (matches_iterator == matches_by_hash.end())
        {
            continue;
        }
        vector <KMerMatch> matches = matches_iterator -> second;
        if (matches.size() > max_k_mer_occurrences)
        {
            continue;
        }
        for (KMerMatch match : matches)
        {
            int relative_shift = i - match.offset;
            size_t overlap = GetSequencesOverlap(seq, sequences[match.seq_number], relative_shift);
            if (overlap < min_overlap)
                continue;
            candidates_set.insert({match.seq_number, relative_shift});
        }
    }

    candidates.resize(0);
    candidates.reserve(candidates_set.size());
    for (auto candidate : candidates_set)
    {
        candidates.push_back(candidate);
    }
}
