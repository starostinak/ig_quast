#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1

#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/modifier.h>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <boost/format.hpp>
#include <mutex>
#include <thread>
#include <chrono>
#include <algorithm>


using namespace seqan;
using std::vector;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::map;
using std::make_pair;
using std::make_tuple;
using bformat = boost::format;

#include "banded_half_smith_waterman.hpp"


template<typename T>
void output_matrix(const T &c) {
  for (const auto &_ : c) {
    for (const auto &__ : _) {
      cout << bformat("%3d ") % __;
    }
    cout << endl;
  }
}


template<typename Ts, typename Tf>
int half_sw_naive(const Ts &s1, const Ts &s2, int indel, int match, int mismatch, const Tf &lizard_tail) {
  vector<vector<int>> c;

  const int INF = 1005000;

  c.resize(length(s1) + 1);
  for (auto &_ : c) {
    _.resize(length(s2) + 1, -INF);
  }

  for (int i1 = 0; i1 <= length(s1); ++i1) {
    c[i1][length(s2)] = lizard_tail(length(s1) - i1);
  }

  for (int i2 = 0; i2 <= length(s2); ++i2) {
    c[length(s1)][i2] = lizard_tail(length(s2) - i2);
  }

  for (int i1 = length(s1) - 1; i1 >= 0; --i1) {
    for (int i2 = length(s2) - 1; i2 >= 0; --i2) {
      int m = (s1[i1] == s2[i2]) ? match : mismatch;
      c[i1][i2] = std::max( { m + c[i1+1][i2+1], indel + c[i1+1][i2], indel + c[i1][i2+1] } );
    }
  }

  output_matrix(c);

  return c[0][0];
}


int main(int argc, char **argv) {
  if (argc < 2) {
    return -1;
  }

  Dna5String s1 = argv[1], s2 = argv[2];

  auto lizard_tail = [](int l) -> int { return ((l > 0) ? -1 : 0); };
  cout << half_sw_naive(s1, s2, -1, 1, -1, lizard_tail) << std::endl;
  cout << half_sw_naive(s2, s1, -1, 1, -1, lizard_tail) << std::endl;
  cout << half_sw_banded(s1, s2, -1, 1, -1, lizard_tail) << std::endl;
  cout << half_sw_banded(s2, s1, -1, 1, -1, lizard_tail) << std::endl;
  cout << half_sw_banded(s2, s1, -1, 1, -1, lizard_tail, 12) << std::endl;
  cout << half_sw_banded(s2, s1, -1, 1, -1, lizard_tail, 18) << std::endl;

  cout << half_sw_banded("MAMA", "AMAM", -1, 1, -1, lizard_tail, 18) << std::endl;
}

