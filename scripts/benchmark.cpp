#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <set>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>

#define LOG_OFFSET		0.01
#define EPSILON			1e-300

using namespace std;

map<string, int> truth_gene, input_gene, truth_barcode, input_barcode;
map<int, string> truth_gene_cvt, input_gene_cvt, truth_barcode_cvt, input_barcode_cvt;
map<int, map<int, double>> truth_matrix, input_matrix;
vector<double> arr_pearson, arr_log_pearson, arr_spearman, arr_emard;
int n_input, n_share, n_truth;
                                                                                                                              
void describe_distribution(vector<double> arr)                                   
{                                                                                
	double average = accumulate(arr.begin(), arr.end(), 0.0)/arr.size();         

	int mid = arr.size() / 2;                             
	int q1 = arr.size() / 4;
	int q3 = arr.size() * 3 / 4;
	sort(arr.begin(), arr.end());
	printf("%.6lf\t%.6lf\t%.6lf\t%.6lf\n", average, arr[mid], arr[q1], arr[q3]); 
}                                                                                

void output_plot()                                                               
{
	cerr << "Score       \tMean     \tMedian  \tQ1      \tQ2      \n";
	cerr << "Spearman    \t";
	describe_distribution(arr_spearman);
	
	cerr << "EMard:      \t";
	describe_distribution(arr_emard);
	cerr << "Pearson:    \t";
	describe_distribution(arr_pearson);
	
	cerr << "Log pearson:\t";
	describe_distribution(arr_log_pearson);
} 

void parse_barcode(string file_path, map<string, int> &barcode_map, map<int, string> &barcode_cvt)
{
	ifstream fi(file_path);
	if (!fi.is_open()) {
		string err = "Could not open file " + file_path;
		perror(err.c_str());
		exit(EXIT_FAILURE);
	}

	string s;
	int id = 0;
	while (fi >> s) {
		++id;
		barcode_map[s.substr(0, 16)] = id;
		barcode_cvt[id] = s.substr(0, 16);
	}

	fi.close();
}

void parse_gene(string file_path, map<string, int> &gene_map, map<int, string> &gene_cvt)
{
	ifstream fi(file_path);
	if (!fi.is_open()) {
		string err = "Could not open file " + file_path;
		perror(err.c_str());
		exit(EXIT_FAILURE);
	}

	string gene_id, gene_name;
    string s;
	int id = 0;
    while (getline(fi, s)) {
        istringstream iss(s);
        iss >> gene_id;
        ++id;
        gene_map[gene_id] = id;
        gene_cvt[id] = gene_id;
    }
	//while (fi >> gene_id >> gene_name) {
//		++id;
//		gene_map[gene_id] = id;
//		gene_cvt[id] = gene_id;
//	}

	fi.close();
}

void parse_matrix(string file_path, map<int, map<int, double>> &matrix_map)
{
	ifstream fi(file_path);
	if (!fi.is_open()) {
		string err = "Could not open file " + file_path;
		perror(err.c_str());
		exit(EXIT_FAILURE);
	}

	string s;
	getline(fi, s);
	while (s[0] == '%') 
		getline(fi, s);

	int gene_id, barcode_id;
	double n_umi;
	do {
		istringstream iss(s);
		iss >> gene_id >> barcode_id >> n_umi;
		matrix_map[barcode_id][gene_id] = n_umi;
	} while(getline(fi, s));

	fi.close();
}

void get_rank(vector<double> &rank, vector<double> &tpm)
{
	vector<pair<double, int>> temp;
	temp.resize(tpm.size());
	int i, j, n = (int)tpm.size();
	
	for (i = 0; i < n; i++) {
		temp[i] = { tpm[i], i };
	}
	sort(temp.begin(), temp.end());
	rank.resize(temp.size());

	for (i = 0; i < n; ) {
		double val, r;
		
		val = temp[i].first;
		j = i + 1;
		while (j < n && val == temp[j].first) {
			j++;
		}

		r = (i + j - 1) / 2.0;
		for ( ; i < j; i++) {
			rank[temp[i].second] = r;
		}
	}
}

double pearson(vector<double> &truth, vector<double> &input)
{
	double mean1 = 0, mean2 = 0;
	double var1 = 0, var2 = 0, num = 0;
	int i = 0, n = truth.size();

	for (i = 0; i < n; i++) {
		mean1 += truth[i];
		mean2 += input[i];
	}
	
	mean1 /= n;
	mean2 /= n;

	for (i = 0; i < n; i++) {
		num += (truth[i] - mean1) * (input[i] - mean2);
		var1 += (truth[i] - mean1) * (truth[i] - mean1);
		var2 += (input[i] - mean2) * (input[i] - mean2);
	}

	if (num == 0 || sqrt(var1 * var2) == 0)
		return 0;
	return num / sqrt(var1 * var2);
}

double log_pearson(vector<double> &truth, vector<double> &input)
{
	double mean1 = 0, mean2 = 0;
	double var1 = 0, var2 = 0, num = 0;
	int i = 0, n = truth.size();

	for (i = 0; i < n; i++) {
		mean1 += log(truth[i] + LOG_OFFSET);
		mean2 += log(input[i] + LOG_OFFSET);
	}
	
	mean1 /= n;
	mean2 /= n;

	for (i = 0; i < n; i++) {
		double x = log(truth[i] + LOG_OFFSET) - mean1;
		double y = log(input[i] + LOG_OFFSET) - mean2;
		num += x * y;
		var1 += x * x;
		var2 += y * y;
	}

	if (num == 0 || sqrt(var1 * var2) == 0)
		return 0;
	return num / sqrt(var1 * var2);
}

double eMARD(vector<double> &truth, vector<double> &input)
{
	double mean = 0;
	int i, n = truth.size();
	for (i = 0; i < n; ++i) {
		if (!(truth[i] + input[i] > EPSILON)) {
			cerr << truth[i] << " " << input[i] << endl;
		}
		assert(truth[i] + input[i] > EPSILON);
		mean += fabs(truth[i] - input[i]) / (truth[i] + input[i]);
	}
	mean /= n;
	return mean;
}

void local_scoring(map<int, double> &truth_cell, map<int, double> &input_cell)
{
	// get all genes express
	map<string, pair<double, double>> temp;
	for (auto &item : truth_cell) {
		if (item.first == 0) {
			assert(item.second == 0);
			continue;
		}
		auto &truth_gene_name = truth_gene_cvt[item.first];
		auto &truth_umi_count = item.second;
		auto &input_umi_count = input_cell[input_gene[truth_gene_name]];
                if (max(truth_umi_count, input_umi_count) <= 1)
                     continue;
		temp[truth_gene_name] = { truth_umi_count, input_umi_count };
		if (truth_umi_count == 0 && input_umi_count == 0) {
			cerr << endl << "1 " << item.first << endl;
		}
		assert(truth_umi_count > 0 || input_umi_count > 0);
	}

	for (auto &item : input_cell) {
		if (item.first == 0) {
			assert(item.second == 0);
			continue;
		}
		auto &input_gene_name = input_gene_cvt[item.first];
		auto &input_umi_count = item.second;
		if (temp.count(input_gene_name)) {
			continue;
		}
		auto &truth_umi_count = truth_cell[truth_gene[input_gene_name]];
                if (max(truth_umi_count, input_umi_count) <= 1)
                     continue;
		temp[input_gene_name] = { truth_umi_count, input_umi_count };
		if (truth_umi_count == 0 && input_umi_count == 0) {
			cerr << endl << "2 " << item.first << endl;
		}
		assert(truth_umi_count > 0 || input_umi_count > 0);
	}

	// generate vector
	vector<double> truth_arr, input_arr;
	for (auto &item : temp) {
		truth_arr.push_back(item.second.first);
		input_arr.push_back(item.second.second);
	}

	// calculate correlation
	vector<double> t_rank, i_rank;
	get_rank(t_rank, truth_arr);
	get_rank(i_rank, input_arr);
	arr_spearman.push_back(pearson(t_rank, i_rank));
	arr_pearson.push_back(pearson(truth_arr, input_arr));
	arr_log_pearson.push_back(log_pearson(truth_arr, input_arr));
	arr_emard.push_back(eMARD(truth_arr, input_arr));
}

void scoring()
{
	for (auto &it : truth_matrix) {
		auto &truth_barcode_id = it.first;
		auto &truth_cell = it.second;
		int input_barcode_id = input_barcode[truth_barcode_cvt[truth_barcode_id]];
		if (input_barcode_id == 0) {
			++n_truth;
			continue;
		}
		++n_share;

		auto &input_cell = input_matrix[input_barcode_id];
		local_scoring(truth_cell, input_cell);
		input_matrix.erase(input_barcode_id);
	}

	n_input = input_matrix.size();
	cerr << "No. shared barcodes: " << n_share << endl;
	cerr << "No. barcodes only in input: " << n_input << endl;
	cerr << "No. barcodes only in truth: " << n_truth << endl;
}

void print_usage()
{
	cerr << "\n";
	cerr << "Usage: ./benchmark_single_cell truth_folder input_folder\n";
	cerr << "Truth folder and input folder must has 3 files: matrix.mtx, genes.tsv, barcodes.tsv\n";	
	cerr << "\n";
}

int main(int argc, char **argv)
{
	if (argc != 3) {
		print_usage();
		return 1;
	}

	parse_gene(string(argv[1]) + "/genes.tsv", truth_gene, truth_gene_cvt);
	parse_gene(string(argv[2]) + "/genes.tsv", input_gene, input_gene_cvt);
	parse_barcode(string(argv[1]) + "/barcodes.tsv", truth_barcode, truth_barcode_cvt);
	parse_barcode(string(argv[2]) + "/barcodes.tsv", input_barcode, input_barcode_cvt);
	parse_matrix(string(argv[1]) + "/matrix.mtx", truth_matrix);
	parse_matrix(string(argv[2]) + "/matrix.mtx", input_matrix);
	scoring();
	output_plot();
	return 0;
}
