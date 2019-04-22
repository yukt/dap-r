#include <regex>
#include <fstream>
#include <cstdio>
#include <string>

using namespace std;

class result_parser
{
public:
	vector<int> model_rank, model_size;
	vector<float> model_posterior, model_score;
	vector<string> model_configure;

	float model_size_mean, model_size_sd;
	float LogNC, Log10NC;
	float min_pip;
	int N;


	vector<int> snp_rank, snp_cluster;
	vector<string> snp_name;
	vector<float> snp_pip, snp_score;

	vector<int> cluster, member;
	vector<float> cluster_pip, cluster_r2;

	result_parser(const char* out_file);
};

result_parser::result_parser(const char* out_file)
{
	regex re_pattern1(" +([^ ]+) +([^ ]+) +([^ ]+) +([^ ]+) +(\\[.+\\]).*");
	string line;
	smatch match;

	ifstream outfile(out_file);

	while (getline(outfile, line))
	{
		regex_search(line, match, re_pattern1);
		if(match.size() != 6) break;
		model_rank.push_back(stoi(match.str(1)));
		model_posterior.push_back(stof(match.str(2)));
		model_size.push_back(stoi(match.str(3)));
		model_score.push_back(stof(match.str(4)));
		model_configure.push_back(match.str(5));
	}

	vector<float> summary_stat; 
	regex re_pattern2("^[A-Z][a-zA-Z =:]*(\\d+\\.*\\d+e?-?\\d*) \\( *[[:alnum:]]+ = ([\\d\\.]+)");
	while (getline(outfile, line))
	{
		regex_search(line, match, re_pattern2);
		if(match.size() != 3) break;
		summary_stat.push_back(stof(match.str(1)));
		summary_stat.push_back(stof(match.str(2)));
	}
	assert(summary_stat.size()==6);
	model_size_mean = summary_stat[0];
	model_size_sd   = summary_stat[1];
	LogNC = summary_stat[2]; 
	Log10NC = summary_stat[3];
	min_pip = summary_stat[4];
	N = (int)summary_stat[5];

	getline(outfile, line);
	getline(outfile, line);

	
	regex re_pattern3("\\(\\((\\d+)\\)\\)[[:blank:]]+([[:graph:]]+)[[:blank:]]+([[:graph:]]+)[[:blank:]]+([[:graph:]]+)[[:blank:]]+([[:graph:]]+)");
	while (getline(outfile, line))
	{
		regex_search(line, match, re_pattern3);
		if(match.size() != 6) break;
		snp_rank.push_back(stoi(match.str(1)));
		snp_name.push_back(match.str(2));
		snp_pip.push_back(stof(match.str(3)));
		snp_score.push_back(stof(match.str(4)));
		snp_cluster.push_back(stoi(match.str(5)));
	}

	getline(outfile, line);
	getline(outfile, line);
	getline(outfile, line);

	
	regex re_pattern4("\\{(\\d+)\\}[[:blank:]]+(\\d+)[[:blank:]]+([[:graph:]]+)[[:blank:]]+([[:graph:]]+)");
	while (getline(outfile, line))
	{
		regex_search(line, match, re_pattern4);
		if(match.size() != 5) break;
		cluster.push_back(stoi(match.str(1)));
		member.push_back(stoi(match.str(2)));
		cluster_pip.push_back(stof(match.str(3)));
		cluster_r2.push_back(stof(match.str(4)));
	}
	outfile.close();
}

