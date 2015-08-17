struct Repertoire {
    struct Cluster {
        size_t id;
        string seq;
    };

    size_t id;
    string filename;
    vector <Cluster> clusters;
    vector <string> sequences;
};

struct NeighbourClusters {
    vector <size_t> ids;
    int score;
};

class RepertoiresWithNeighbours {
public:
    RepertoiresWithNeighbours(const Repertoire & rep1, const Repertoire & rep2,
                              const unordered_map <size_t, NeighbourClusters> & neighbours12,
                              const unordered_map <size_t, NeighbourClusters> & neighbours21)
                          : rep1_(rep1), rep2_(rep2), neighbours12_(neighbours12), neighbours21_(neighbours21) {}

    void Write(ostream & out) const;
private:
    static void WriteNeighbours(ostream & out, size_t rep_id1, size_t rep_id2,
                                 unordered_map <size_t, NeighbourClusters> const & neighbours);

    const Repertoire & rep1_;
    const Repertoire & rep2_;
    unordered_map <size_t, NeighbourClusters> neighbours12_;
    unordered_map <size_t, NeighbourClusters> neighbours21_;
};

void RepertoiresWithNeighbours::Write(ostream & out) const
{
    WriteNeighbours(out, rep1_.id, rep2_.id, neighbours12_);
    WriteNeighbours(out, rep2_.id, rep1_.id, neighbours21_);
}

void RepertoiresWithNeighbours::WriteNeighbours(ostream & out, size_t rep_id1, size_t rep_id2,
                                            const unordered_map<size_t, NeighbourClusters> & neighbours) {
    map <size_t, NeighbourClusters> sorted(neighbours.begin(), neighbours.end());
    for (auto it = sorted.begin(); it != sorted.end(); ++it) {
        size_t first_cid = it->first;
        int score = it->second.score;
        for (auto second_cid: it->second.ids) {
            out << rep_id1 << '.' << first_cid << " - " << rep_id2 << '.' << second_cid <<
                " (score = " << score << ")" << endl;
        }
    }
}
