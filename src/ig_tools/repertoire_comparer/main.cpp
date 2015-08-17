#include "../utils/fasta_reader.hpp"
#include "../utils/sequence_tools.hpp"
#include "../utils/string_tools.hpp"
#include "neighbours_constructor.hpp"
#include "getopt.h"

#define HASH_SIZE 80
#define INDENT 150
#define MATCH_CUTOFF 2000
#define MIN_MATCH_CUTOFF 100

struct CmdParams {
    int threads_num;
    vector <string> input_files;
    string output_file;
    int max_mismatches;
    int max_gaps;
    int min_overlap;
    int match_score;
    int mismatch_penalty;
    int gap_penalty;
};

int parse_cmd_parameters(int argc, char *argv[], CmdParams &params) {
    const char *short_options = "";
    const struct option long_options[] = {
        {"out", required_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'},
        {"max_mismatches", required_argument, NULL, 'm'},
        {"max_gaps", required_argument, NULL, 'g'},
        {"min_overlap", required_argument, NULL, 'l'},
        {"match_score", required_argument, NULL, 's'},
        {"mismatch_penalty", required_argument, NULL, 'i'},
        {"gap_penalty", required_argument, NULL, 'p'},
        {NULL, 0, NULL, 0}
    };
    int rez = 0;
    int option_index = 0;
    params.match_score = 1;
    params.gap_penalty = 10;
    params.mismatch_penalty = 5;
    params.max_gaps = 0;
    params.max_mismatches = 4;
    params.min_overlap = 300;
    params.threads_num = 1;
    while ((rez = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1) {
        switch(rez) {
        case 'o':
            params.output_file = optarg;
            break;
        case 't':
            params.threads_num = atoi(optarg);
            break;
        case 'm':
            params.max_mismatches = atoi(optarg);
            break;
        case 'g':
            params.max_gaps = atoi(optarg);
            break;
        case 'l':
            params.min_overlap = atoi(optarg);
            break;
        case 's':
            params.match_score = atoi(optarg);
            break;
        case 'i':
            params.mismatch_penalty = atoi(optarg);
            break;
        case 'p':
            params.gap_penalty = atoi(optarg);
            break;
        default:
            cerr << "Error: unknown parameter - " << rez << endl;
            return -1;
        }
    }
    for (int i = optind; i != argc; ++i) {
        params.input_files.push_back(argv[i]);
    }
    if (params.output_file.empty() || params.input_files.size() < 2) {
        cerr << "You must specify output file and at least two input files" << endl;
        return -1;
    }
    return 0;
}

int read_repertoires(vector <string> const & filenames, vector <Repertoire> & repertoires) {
    for (size_t i = 0; i != filenames.size(); ++i) {
        FastaReader <fasta_read, FastaReadConstructor> reader(filenames[i]);
        vector <fasta_read> reads = reader.Read();
        repertoires[i].id = i + 1;
        repertoires[i].filename = filenames[i];
        for (auto read_it = reads.begin(); read_it != reads.end(); ++read_it) {
            vector <string> fields = split(read_it->name, "___");
            if (fields.size() != 4) {
                return -1;
            }
            size_t cluster_id = atoi(fields[1].c_str());
            repertoires[i].clusters.push_back({cluster_id, read_it->seq});
            repertoires[i].sequences.push_back(read_it->seq);
        }
    }
    return 0;
}

void write_neighbour_clusters(ostream & out, vector <Repertoire> const & repertoires, vector <RepertoiresWithNeighbours> const & neighbour_clusters) {
    for (size_t i = 0; i != repertoires.size(); ++i) {
        out << repertoires[i].id << '.' << repertoires[i].filename << ' ';
    }
    out << endl;

    for (auto repsWithNeighbours : neighbour_clusters) {
        repsWithNeighbours.Write(out);
    }
}

void print_usage() {
    cout << "./repertoire_comparer errors outfile in1.fa in2.fa ..." << endl;
}

int main(int argc, char ** argv) {
    // 4 4 out small1.fasta small2.fasta
    CmdParams cmd_params;
    if (parse_cmd_parameters(argc, argv, cmd_params) != 0) {
        return -1;
    }
    vector <Repertoire> repertoires(cmd_params.input_files.size());
    if (read_repertoires(cmd_params.input_files, repertoires) < 0) {
        cerr << "Repertoire file formatted incorrectly" << endl;
        return -1;
    }

    vector <RepertoiresWithNeighbours> neighbour_clusters;
    
    for (size_t i = 0; i != repertoires.size(); ++i) {
        for (size_t j = i + 1; j != repertoires.size(); ++j) {
            int seq_match_cutoff = max(MIN_MATCH_CUTOFF, 
                                       int(min(repertoires[i].clusters.size(), repertoires[j].clusters.size()) * 0.1));
            ScoringParameters scoring_params;
            scoring_params.gap_penalty = cmd_params.gap_penalty;
            scoring_params.match_score = cmd_params.match_score;
            scoring_params.max_gaps = cmd_params.max_gaps;
            scoring_params.max_mismatches = cmd_params.max_mismatches;
            scoring_params.min_overlap = cmd_params.min_overlap;
            scoring_params.mismatch_penalty = cmd_params.mismatch_penalty;
            unique_ptr <SequenceScorer> scorer;
            if (scoring_params.max_gaps == 0) {
                scorer = unique_ptr<SequenceScorer>(new HammingSequenceScorer(scoring_params, HASH_SIZE));
            } else {
                scorer = unique_ptr<SequenceScorer>(new SWSequenceScorer(scoring_params, HASH_SIZE));
            }
            NeighboursConstructor constructor(HASH_SIZE, cmd_params.min_overlap, min(seq_match_cutoff, MATCH_CUTOFF), INDENT,
                                              *scorer, cmd_params.threads_num);
            neighbour_clusters.push_back(constructor.Construct(repertoires[i], repertoires[j]));
        }
    }

    ofstream fout(cmd_params.output_file);
    write_neighbour_clusters(fout, repertoires, neighbour_clusters);
    fout.close();

    return 0;
}
