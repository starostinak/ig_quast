#include "../utils/include_me.hpp"
#include "../utils/fastq_reader.hpp"
#include "../utils/cluster_utils/clusterizator.hpp"
#include "../utils/string_tools.hpp"
#include "omp.h"
#include "../utils/cluster_utils/clusters.hpp"
#include "clusters_evaluator.hpp"

string GetResultFname(string reads_fname) {
	size_t pos = reads_fname.find(".fastq");
	string base_name = reads_fname.substr(0, reads_fname.size() - 6);
	if(pos == string::npos) {
		pos = reads_fname.find(".fq");
		base_name = reads_fname.substr(0, reads_fname.size() - 3);
	}
	return base_name + ".clusters";
}

string CreateBaseName(string fname) {
	string suffix = ".clusters";
	return fname.substr(0, fname.size() - suffix.size());
}

int main(int argc, char *argv[]) {
	/*
	 * argv[1] - ideal clusters FA
	 * argv[2] - ideal clusters RCM
	 * argv[3] - clusters to evaluation FA
	 * argv[4] - clusters to evaluation RCM
	 */

	if(argc < 3) {
		cout << "Invalid output parameters" << endl <<
				"\targv[1] - file with .rcm for original repertoire (mandatory)" << endl <<
				"\targv[2] - file with .clusters for original repertoire (mandatory)" << endl <<
				"\targv[3] - file with .rcm for constructed repertoire (mandatory)" << endl <<
				"\targv[4] - file with .clusters for constructed repertoire (mandatory)" << endl <<
				"\targv[5] - basename of resulting files";
		return 1;
	}

	ClustersEvaluator evaluator(argv[1], argv[2], argv[3], argv[4]);
	evaluator.Evaluate();
	evaluator.GetMetrics().Print(cout);
	string basename = CreateBaseName(string(argv[3]));
	if(argc > 3)
		basename = string(argv[5]);
	string metrics_fname = basename + ".txt";
	ofstream out(metrics_fname.c_str());
	cout << "Metrics were written to " << metrics_fname << endl;
	evaluator.GetMetrics().Print(out);
	out.close();

	// csv
	metrics_fname = basename + ".csv";
	ofstream out2(metrics_fname.c_str());
	evaluator.GetMetrics().PrintCSV(out2);
	out2.close();

	// size/percentage output
	string length_perc_fname = basename + ".size.perc";
	ofstream out3(length_perc_fname.c_str());
	evaluator.PrintSizeIdentityPercentage(out3);
	out3.close();

	//ClusterSizeHistogramFiller(evaluator.OriginalClusters()).WriteToFile(basename + ".original_clusters_sizes.txt");
	//ClusterSizeHistogramFiller(evaluator.ConstructedClusters()).WriteToFile(basename + ".constructed_clusters_sizes.txt");
}
