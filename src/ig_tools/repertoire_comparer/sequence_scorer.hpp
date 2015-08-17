#include "../../include/ssw/ssw_cpp.h"
#include "seqan/align.h"
#include "striped_smith_waterman.hpp"

struct ScoringParameters {
    int max_mismatches;
    int max_gaps;
    int min_overlap;
    int match_score;
    int mismatch_penalty;
    int gap_penalty;
};

class SequenceScorer {
public:
    SequenceScorer(ScoringParameters const &params, int kmer_length);

    virtual int score_sequences(string const & seq1, size_t kmer1_start,
                        string const & seq2, size_t kmer2_start, int curr_best_score) const = 0;

protected:
    ScoringParameters params_;
    int kmer_length_;
};

class SWSequenceScorer: public SequenceScorer {
public:
    SWSequenceScorer(ScoringParameters const &params, int kmer_length)
        : SequenceScorer(params, kmer_length) { }
    virtual int score_sequences(string const & seq1, size_t kmer1_start,
                        string const & seq2, size_t kmer2_start, int curr_best_score) const;
};

class HammingSequenceScorer: public SequenceScorer {
public:
    HammingSequenceScorer(ScoringParameters const &params, int kmer_length)
        : SequenceScorer(params, kmer_length) { }
    virtual int score_sequences(string const & seq1, size_t kmer1_start,
                        string const & seq2, size_t kmer2_start, int curr_best_score) const;
};

SequenceScorer::SequenceScorer(ScoringParameters const &params, int kmer_length)
    : params_(params), kmer_length_(kmer_length) { }

int SWSequenceScorer::score_sequences(const string &seq1, size_t start1, const string &seq2, size_t start2, int curr_best_score) const {
    using namespace seqan;
    typedef String<char> TSequence;
    typedef Align<TSequence, ArrayGaps> TAlign;
    TSequence seq1_cut = seq1.substr(start1);
    TSequence seq2_cut = seq2.substr(start2);

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq1_cut);
    assignSource(row(align, 1), seq2_cut);

    Score <int, Simple> scoring_scheme(params_.match_score, -params_.mismatch_penalty, -params_.gap_penalty);

    int score = globalAlignment(align, scoring_scheme, AlignConfig <true, true, true, true>(), -params_.max_gaps, params_.max_gaps);

    if (score < curr_best_score) {
        return 0;
    }

    AlignmentStats stats;
    computeAlignmentStats(stats, align, scoring_scheme);
    if (stats.numMatches + stats.numMismatches < params_.min_overlap
            || stats.numMismatches > params_.max_mismatches
            || stats.numGapOpens + stats.numGapExtensions > params_.max_gaps) {
        return 0;
    }

    return score;
}


int HammingSequenceScorer::score_sequences(const std::string &seq1, size_t kmer1_start, const std::string &seq2, size_t kmer2_start, int curr_best_score) const {
    size_t seq1_start = (kmer1_start < kmer2_start) ? 0 : kmer1_start - kmer2_start;
    size_t seq2_start = (kmer1_start < kmer2_start) ? kmer2_start - kmer1_start : 0;
    size_t overlap_length = min(seq1.size() - seq1_start, seq1.size() - seq2_start);
    int max_possible_score = overlap_length * params_.match_score;
    if (overlap_length * params_.match_score < curr_best_score) {
        return 0;
    }
    int score = 0;
    for (size_t i = 0; i != overlap_length; ++i) {
        if (seq1[seq1_start + i] == seq2[seq2_start + i]) {
            score += params_.match_score;
        } else {
            score -= params_.mismatch_penalty;
            max_possible_score -= params_.mismatch_penalty;
        }
        if (max_possible_score < curr_best_score) {
            return 0;
        }
    }
    if (score < curr_best_score) {
        return 0;
    }
    return score;
}
