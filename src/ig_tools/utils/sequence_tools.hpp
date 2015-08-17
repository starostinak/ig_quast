#pragma once

#include "include_me.hpp"

char get_complementary(char nucl){
	if(nucl == 'A' || nucl == 'a')
		return 'T';
	if(nucl == 'T' || nucl == 't')
		return 'A';
	if(nucl == 'C' || nucl == 'c')
		return 'G';
	if(nucl == 'G' || nucl == 'g')
		return 'C';
	if(nucl == 'N' || nucl == 'n')
		return 'A';
	cout << "Char " << nucl << " is not from nucleotide alphabet" << endl;
	assert(false);
	return 'A';
}

string reverse_complementary(string seq) {
	string rc_seq = seq;
	for(size_t i = 0; i < seq.size(); i++)
		rc_seq[seq.size() - i - 1] = get_complementary(seq[i]);
	return rc_seq;
}

char get_complementary_nn(char nucl){
	if(nucl == 'A' || nucl == 'a')
		return 'T';
	if(nucl == 'T' || nucl == 't')
		return 'A';
	if(nucl == 'C' || nucl == 'c')
		return 'G';
	if(nucl == 'G' || nucl == 'g')
		return 'C';
	if(nucl == 'N' || nucl == 'n')
		return 'N';
	cout << "Char " << nucl << " is not from nucleotide alphabet" << endl;
	assert(false);
	return 'N';
}

string reverse_complementary_nn(string seq) {
	string rc_seq(seq.size(), 'N');
	transform(seq.begin(), seq.end(), rc_seq.begin(), get_complementary_nn);
	reverse(rc_seq.begin(), rc_seq.end());
	return rc_seq;
}

size_t HammingDistance(string s1, string s2) {
	assert(s1.size() == s2.size());
	size_t dist = 0;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			dist++;
	return dist;
}

set<size_t> DifferentPositions(string s1, string s2) {
	assert(s1.size() == s2.size());
	set<size_t> pos;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			pos.insert(i);
	return pos;
}

string random_correction(string s1, string s2) {
	assert(s1.size() == s2.size());
	string s = s1;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			s[i] = s2[i];
	return s;
}

char GetRangomNucleotide() {
	static char nucls[4] = { 'A', 'C', 'G', 'T'};
	return nucls[rand() % 4];
}

string GetPalindrom(size_t half_size) {
	string str;
	for(size_t i = 0; i < half_size; i++)
		str = str + GetRangomNucleotide();
	for(size_t i = 0; i < half_size; i++)
		str = str + get_complementary(str[half_size - i - 1]);
	return str;
}

char GetAnotherRandomNucleotide(char nucl) {
	char res;
	do
		res = GetRangomNucleotide();
	while(res == nucl);
	return res;
}

int get_string_mismatches(string const & seq1, size_t start1,
								  string const & seq2, size_t start2,
								  int max_errors) {
	int errors = 0;
	size_t size1 = seq1.size();
	size_t size2 = seq2.size();
	size_t len = min(size1 - start1, size2 - start2);
	if (max_errors == 0) {
		if (seq1.compare(start1, len, seq2, start2, len) == 0) {
			return 0;
		}
		return len;
	}
	// for (size_t i = 0; i + start1 != size1 && i + start2 != size2; ++i) {
	char const * c1 = seq1.c_str() + start1;
	char const * c2 = seq2.c_str() + start2;
	for (size_t i = 0; i != len; ++i) {
		//if (seq1[i + start1] != seq2[i + start2]) {
		if (*c1 != *c2) {
			++errors;
		}
		if (errors > max_errors) {
			return len;
		}
		++c1; ++c2;
	}
	return errors;
}

size_t score_sequences(string const & seq1, size_t start1, string const & seq2, size_t start2, int max_errors, int error_penalty) {
	size_t len = min(seq1.size() - start1, seq2.size() - start2);
	size_t errors = get_string_mismatches(seq1, start1, seq2, start2, max_errors);
	if (errors > max_errors) {
		return 0;
	}
	if (len < errors * error_penalty) {
		return 0;
	}
	return len - get_string_mismatches(seq1, start1, seq2, start2, max_errors) * error_penalty; 
}

size_t get_string_matches(string const & seq1, string const & seq2, int max_errors, int max_shift) {
	int errors = max_errors;
	int prev_errors = errors;
	size_t matches = min(seq1.size(), seq2.size()) - errors;
	for (int shift = 0; shift != max_shift; ++shift) {
		size_t best_matches = max(min(seq1.size(), seq2.size() - shift), min(seq1.size() - shift, seq2.size()));
		if (best_matches <= matches) {
			break;
		}
		errors = min(get_string_mismatches(seq1, 0, seq2, shift, errors), errors);
		errors = min(get_string_mismatches(seq1, 0, reverse_complementary_nn(seq2), shift, errors), errors);
		if (errors != prev_errors) {
			matches = max(matches, min(seq1.size(), seq2.size() - shift) - errors);
		}
		if (shift != 0) {
			errors = min(get_string_mismatches(seq1, shift, seq2, 0, errors), errors);
			errors = min(get_string_mismatches(seq1, shift, reverse_complementary_nn(seq2), 0, errors), errors);
		}
		if (errors != prev_errors) {
			matches = max(matches, min(seq1.size() - shift, seq2.size()) - errors);
		}
	}
	return matches;
}
