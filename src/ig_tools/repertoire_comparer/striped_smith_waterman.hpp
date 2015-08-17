#include <vector>
#include <string>
#include <cmath>
using std::vector;
using std::string;
using std::max;

class Aligner {
public:
    Aligner(int match_score, int mismatch_penalty, int gap_penalty, int max_mismatches, int max_gaps);

    int Align(string const &s1, string const &s2, int best_known_score) const;

private:
    struct Cell {
        Cell(): valid(false), score(0), mismatches(0), gaps(0) { }

        bool valid;
        int score;
        int mismatches;
        int gaps;
    };

    void convert_matrix_coords(int &i, int &j) const;

    Cell get_cell(vector <vector <Cell> > const &matrix, int i, int j) const;

    int match_score_;
    int mismatch_penalty_;
    int gap_penalty_;
    int max_mismatches_;
    int max_gaps_;
};


Aligner::Aligner(int match_score, int mismatch_penalty, int gap_penalty, int max_mismatches, int max_gaps)
    : match_score_(match_score), mismatch_penalty_(mismatch_penalty), gap_penalty_(gap_penalty), max_mismatches_(max_mismatches), max_gaps_(max_gaps) {
}

void Aligner::convert_matrix_coords(int &i, int &j) const {
    i = i - j + max_gaps_;
}

Aligner::Cell Aligner::get_cell(const vector<vector<Aligner::Cell> > &matrix, int i, int j) const {
    convert_matrix_coords(i, j);
    if (i < 0 || j < 0 || i >= matrix.size()) {
        return Cell();
    }
    return matrix[i][j];
}


int Aligner::Align(const string &s1, const string &s2, int best_known_score) const {
    vector <vector <Cell> > matrix(max_gaps_ * 2 + 1, vector <Cell>(s1.size() + 1));
    matrix[max_gaps_][0].valid = true;
    for (size_t i = 1; i != max_gaps_; ++i) {
        matrix[i + max_gaps_][0].score = matrix[i + max_gaps_ - 1][0].score + gap_penalty_;
        matrix[i + max_gaps_][0].gaps = matrix[i + max_gaps_ - 1][0].gaps + 1;
        matrix[i + max_gaps_][0].valid = true;
        matrix[max_gaps_ - i][i].score = matrix[max_gaps_ - i + 1][i - 1].score + gap_penalty_;
        matrix[max_gaps_ - i][i].gaps = matrix[max_gaps_ - i + 1][i - 1].gaps + 1;
        matrix[max_gaps_ - i][i].valid = true;
    }

    for (int i1 = 0; i1 != static_cast<int>(s1.size()); ++i1) {
        bool any_valid = false;
        int best_score = 0;
        for (int i2 = i1 - max_gaps_; i2 <= i1 + max_gaps_ && i2 != static_cast<int>(s2.size()); ++i2) {
            int i = i2 + 1;
            int j = i1 + 1;
            Cell diag = get_cell(matrix, i - 1, j - 1);
            Cell up = get_cell(matrix, i, j - 1);
            Cell left = get_cell(matrix, i - 1, j);

            convert_matrix_coords(i, j);

            if (diag.valid) {
                matrix[i][j].valid = true;
                matrix[i][j].score = diag.score;
                matrix[i][j].mismatches = diag.mismatches;
                if (s1[i1] != s2[i2]) {
                    matrix[i][j].score -= mismatch_penalty_;
                    matrix[i][j].mismatches += 1;
                    if (matrix[i][j].mismatches > max_mismatches_) {
                        matrix[i][j].valid = false;
                    }
                } else {
                    matrix[i][j].score += match_score_;
                }
            }

            if (up.valid && up.gaps < max_gaps_ &&
                    (up.score - gap_penalty_ > matrix[i][j].score || !matrix[i][j].valid)) {
                matrix[i][j].valid = true;
                matrix[i][j].score = up.score - gap_penalty_;
                matrix[i][j].gaps = up.gaps + 1;
            }

            if (left.valid && left.gaps < max_gaps_ &&
                    (left.score - gap_penalty_ > matrix[i][j].score || !matrix[i][j].valid)) {
                matrix[i][j].valid = true;
                matrix[i][j].score = left.score - gap_penalty_;
                matrix[i][j].gaps = left.gaps + 1;
            }

            any_valid |= matrix[i][j].valid;
            best_score = max(best_score, matrix[i][j].score);
        }
        if (!any_valid || best_score + match_score_ * (s1.size() - i1 - 1) < best_known_score) {
            break;
        }
    }
    int best_score = 0;
    for (size_t i = 0; i != max_gaps_ + 1; ++i) {
        if (matrix[i][s1.size()].valid && matrix[i][s1.size()].score > best_score) {
            best_score = matrix[i][s1.size()].score;
        }
    }
    for (size_t i = 1; i != max_gaps_ + 1; ++i) {
        if (matrix[max_gaps_ + i][s1.size() - i].valid && matrix[max_gaps_ + i][s1.size() - i].score > best_score) {
            best_score = matrix[max_gaps_ + i][s1.size() - i].score;
        }
    }
    if (best_score > best_known_score) {
        return best_score;
    }
    return 0;
}
