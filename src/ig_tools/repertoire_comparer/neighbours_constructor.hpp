#pragma once

#include <unordered_set>
#include <unordered_map>
#include <functional>

#include <boost/thread/thread.hpp>
#include <boost/asio/io_service.hpp>

#include "../utils/sequence_tools.hpp"
#include "shared_kmer_string_set.hpp"
#include "sequence_scorer.hpp"
#include "repertoires.hpp"

class FirstAndLastIndexSelector : public IndexSelector {
public:
    FirstAndLastIndexSelector(size_t indent) : indent_(indent) {}

    virtual bool NeedIncludeToIndex(const string & seq, size_t offset, size_t length) const
    {
        return offset < indent_ || offset + length > seq.length() - indent_;
    }

private:
    const size_t indent_;
};

class NeighboursConstructor {
public:
    NeighboursConstructor(size_t hash_size, int min_overlap, size_t match_cutoff, size_t indent,
                          const SequenceScorer &seq_scorer, size_t threads_num);

    RepertoiresWithNeighbours Construct(const Repertoire & rep1, const Repertoire & rep2) const;

private:
    struct Adjacency {
        size_t id;
        string const * seq;
        bool strand;
        int shift;
    };

    void init_neighbours(const vector <Repertoire::Cluster> & clusters, unordered_map <size_t, NeighbourClusters> & neighbours) const;
    void update_NeighbourClusters(const Repertoire::Cluster & cluster,
                                   const Repertoire & repertoire,
                                   const SharedKMerStringSet & indexedRepertoire,
                                   unordered_map <size_t, NeighbourClusters> & neighbours12,
                                   unordered_map <size_t, NeighbourClusters> & neighbours21) const;
    void update_score(Repertoire::Cluster const & cluster, Adjacency const & adj, NeighbourClusters & neighbour12_seq, NeighbourClusters & neighbour21_seq) const;
    static size_t GetSequencesOverlap(const string & seq1, const string & seq2, int relative_shift);

    size_t hash_size_;
    int min_overlap_;
    size_t match_cutoff_;
    size_t indent_;
    size_t threads_num_;
    SequenceScorer const & seq_scorer_;
    mutable boost::shared_mutex mutex_;
};

NeighboursConstructor::NeighboursConstructor(size_t hash_size, int min_overlap, size_t match_cutoff, size_t indent,
                                             const SequenceScorer & seq_scorer, size_t threads_num)
        : hash_size_(hash_size), min_overlap_(min_overlap), match_cutoff_(match_cutoff), indent_(indent),
          seq_scorer_(seq_scorer), threads_num_(threads_num) { }

RepertoiresWithNeighbours NeighboursConstructor::Construct(const Repertoire & rep1, const Repertoire & rep2) const {
    auto indexedRep2 = SharedKMerStringSet(std::shared_ptr<IndexSelector>(new FirstAndLastIndexSelector(indent_)),
                                           hash_size_,
                                           min_overlap_,
                                           match_cutoff_);
    for (auto seq_it = rep2.sequences.begin(); seq_it != rep2.sequences.end(); ++seq_it) {
        indexedRep2.Insert(*seq_it);
    }
    size_t max_cluster_id = 0;
    for (auto rep2Cluster : rep2.clusters) {
        max_cluster_id = max(rep2Cluster.id, max_cluster_id);
    }

    /*
    cout << ".. Processing " << rep1.filename << " - " << rep2.filename << endl;
    cout << ".. Processed 0%" << std::flush;
    */

    unordered_map <size_t, NeighbourClusters> neighbours12;
    unordered_map <size_t, NeighbourClusters> neighbours21;
    init_neighbours(rep1.clusters, neighbours12);
    init_neighbours(rep2.clusters, neighbours21);

    boost::asio::io_service io_service;
    boost::thread_group threadpool;
    {
        boost::asio::io_service::work work(io_service);
        for (size_t i = 0; i != threads_num_; ++i) {
            threadpool.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
        }

        for (auto it = rep1.clusters.begin(); it != rep1.clusters.end(); ++it) {
            io_service.post(boost::bind(&NeighboursConstructor::update_NeighbourClusters,
                                        this, cref(*it), cref(rep2), cref(indexedRep2),
                                        ref(neighbours12), ref(neighbours21)));
        }
    }
    threadpool.join_all();

    return RepertoiresWithNeighbours(rep1, rep2, neighbours12, neighbours21);
}

void NeighboursConstructor::init_neighbours(const vector <Repertoire::Cluster> & clusters,
                                            unordered_map<size_t, NeighbourClusters> &neighbours) const {
    for (auto cluster : clusters) {
        neighbours[cluster.id] = { {}, 0};
    }
}

void NeighboursConstructor::update_NeighbourClusters(const Repertoire::Cluster & cluster,
                                                      const Repertoire & repertoire,
                                                      const SharedKMerStringSet & indexedRepertoire,
                                                      unordered_map <size_t, NeighbourClusters> & neighbours12,
                                                      unordered_map <size_t, NeighbourClusters> & neighbours21) const {
    NeighbourClusters & neighbour12_seq = neighbours12[cluster.id];
    neighbour12_seq = { {}, 0};
    vector <SharedKMerString> matchCandidates;
    indexedRepertoire.Find(cluster.seq, matchCandidates);
    for (auto candidate : matchCandidates)
    {
        const Repertoire::Cluster & neighbour_cluster = repertoire.clusters[candidate.seq_number_in_set];
        Adjacency adj_seq = {neighbour_cluster.id, &neighbour_cluster.seq, true, candidate.pattern_start_position};
        // TODO strand always is true
        update_score(cluster, adj_seq, neighbour12_seq, neighbours21[neighbour_cluster.id]);
    }
}

void NeighboursConstructor::update_score(const Repertoire::Cluster & cluster, Adjacency const & adj,
                                         NeighbourClusters & neighbour12_seq, NeighbourClusters & neighbour21_seq) const {
    size_t start1 = (adj.shift < 0) ? 0 : adj.shift;
    size_t start2 = (adj.shift < 0) ? -adj.shift : 0;
    /*
    size_t max_score = GetSequencesOverlap(cluster.seq, *adj.seq, adj.shift);
    size_t min_known_score = min(neighbour12_seq.score, neighbour21_seq.score);
    if (max_score < min_known_score) {
        return;
    }
    int curr_max_errors = (max_score - min_known_score) / error_penalty_;
    size_t score = score_sequences(cluster.seq, start1, *adj.seq, start2,
                                   min(max_errors_, curr_max_errors), error_penalty_);
    */
    size_t overlap_len = GetSequencesOverlap(cluster.seq, *adj.seq, adj.shift);
    if (overlap_len < min_overlap_) {
        return;
    }
    size_t best_known_score = min(neighbour12_seq.score, neighbour21_seq.score);
    int score = seq_scorer_.score_sequences(cluster.seq, start1, *adj.seq, start2, best_known_score);
    if (score == 0) {
        return;
    }

    if (score > neighbour12_seq.score) {
        neighbour12_seq = { { adj.id }, score };
    } else if (score == neighbour12_seq.score) {
        neighbour12_seq.ids.push_back(adj.id);
    }

    {
        boost::upgrade_lock<boost::shared_mutex> lock(mutex_);
        if (score > neighbour21_seq.score) {
            boost::upgrade_to_unique_lock<boost::shared_mutex> unique_lock(lock);
            neighbour21_seq = { { cluster.id }, score };
        } else if (score == neighbour21_seq.score) {
            boost::upgrade_to_unique_lock<boost::shared_mutex> unique_lock(lock);
            neighbour21_seq.ids.push_back(cluster.id);
        }
    }
}

size_t NeighboursConstructor::GetSequencesOverlap(const string & seq1, const string & seq2, int relative_shift)
{ // TODO move this function to external utils
    if (relative_shift < 0) {
        return min(seq1.length(), seq2.length() + relative_shift);
    } else {
        return min(seq1.length() - relative_shift, seq2.length());
    }
}
