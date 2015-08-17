#pragma once

#include "../include_me.hpp"
#include "../string_tools.hpp"
#include "../fasta_reader.hpp"

class Clusters {
	ifstream rcm_src_;
	string fa_filename_;

	map<string, size_t> read_cluster_;
	map<string, size_t> read_index_;
	vector<string> reads_;
	map<size_t, vector<size_t> > clusters_;
	map<size_t, string> cluster_seqs_;

public:
	Clusters(string rcm_fname) :
		rcm_src_(rcm_fname.c_str()), fa_filename_() {
		assert(!rcm_src_.fail());
	}
	Clusters(string rcm_fname, string fa_fname) :
		rcm_src_(rcm_fname.c_str()), fa_filename_(fa_fname) {
		assert(!rcm_src_.fail());
	}

	void ExtractFromFile() {
		while(!rcm_src_.eof()) {
			string tmp;
			getline(rcm_src_, tmp);
			if(tmp != "") {
				vector<string> splits = split(tmp, '\t');
				if(splits.size() != 2) {
					cout << "String " << tmp << " is wrong. Number of splits is " << splits.size() << endl;
					assert(splits.size() == 2);
				}
				string name;
				size_t pos = splits[0].find("fr=");
				name = splits[0].substr(0, pos);
				read_cluster_[/*splits[0]*/name] = string_to_number<size_t>(splits[1]);
			}
		}

		if (!fa_filename_.empty()) {
			vector<fasta_read> seqs = FastaReader<fasta_read, FastaReadConstructor>(fa_filename_).Read();
			for(auto it = seqs.begin(); it != seqs.end(); ++it) {
				vector<string> fields = split(it->name, "___");
				if(fields.size() != 4) {
					cout << "Fasta cluster description " << it->name << " is wrong. Format must be '>cluster___id___size___sz'" << endl;
					assert(fields.size() == 4);
				}
				cluster_seqs_[atol(fields[1].c_str())] = it->seq;
			}
		}

		size_t index = 0;
		for(auto it = read_cluster_.begin(); it != read_cluster_.end(); it++) {
			read_index_[it->first] = index;
			index++;
			reads_.push_back(it->first);
		}

		for(auto it = read_cluster_.begin(); it != read_cluster_.end(); it++) {
			clusters_[it->second].push_back(read_index_[it->first]);
		}
	}

	vector<string> Reads() { return reads_; }

	map<size_t, vector<size_t> >::iterator clusters_begin() { return clusters_.begin(); }

	map<size_t, vector<size_t> >::iterator clusters_end() { return clusters_.end(); }

	size_t GetCluster(size_t read_index) {
		assert(read_index < reads_.size());
		return read_cluster_[reads_[read_index]];
	}

	size_t GetCluster(string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return GetCluster(read_index_[read_name]);
	}

	size_t GetClusterSizeByReadName(string read_name) {
		if(read_index_.find(read_name) == read_index_.end())
			return size_t(-1);
		return ClusterSize(GetCluster(read_name));
	}

	size_t ClusterSize(size_t cluster) {
		return clusters_[cluster].size();
	}

	vector<size_t> GetReadsByCluster(size_t cluster) {
		return clusters_[cluster];
	}

	size_t ClustersNumber() {
		return clusters_.size();
	}

	bool IsClusterSingleton(size_t index) {
		return clusters_[index].size() == 1;
	}

	bool ReadExists(string read_name) {
		return read_index_.find(read_name) != read_index_.end();
	}

	string const & GetClusterSequence(size_t cluster) {
		return cluster_seqs_[cluster];
	}

};
