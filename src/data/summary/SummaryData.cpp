#include "data/summary/SummaryData.h"

#include "infrastructure/Logger.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <zlib.h>

using std::map;
using std::string;
using std::vector;

namespace {

void for_each_line(const string &filename, const std::function<void(int, const string &)> &callback) {
    int line_number = 0;
    if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
        gzFile input = gzopen(filename.c_str(), "rb");
        if (!input) LOGGER.e(0, "cannot open the file [" + filename + "] to read.");
        char buffer[1024];
        while (gzgets(input, buffer, sizeof(buffer))) {
            string line(buffer);
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            callback(++line_number, line);
        }
        gzclose(input);
        return;
    }
    std::ifstream input(filename);
    if (!input) LOGGER.e(0, "cannot open the file [" + filename + "] to read.");
    string line;
    while (std::getline(input, line)) callback(++line_number, line);
}

vector<string> fields(const string &line) {
    std::stringstream stream(line);
    vector<string> values;
    string value;
    while (stream >> value) values.push_back(value);
    return values;
}

string uppercase(string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::toupper(c); });
    return value;
}

double median(vector<double> values) {
    if (values.empty()) return 0.0;
    std::stable_sort(values.begin(), values.end());
    size_t size = values.size();
    if (size % 2) return values[(size - 1) / 2];
    return (values[size / 2] + values[size / 2 - 1]) / 2.0;
}

vector<string> read_snp_metafile(const string &filename, map<string, int> &gws_snp_name_map, double thresh) {
    vector<string> snplist;
    for_each_line(filename, [&](int line_number, const string &line) {
        vector<string> values = fields(line);
        if (values.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + filename + "] should be in GCTA-COJO format, line " + std::to_string(line_number) + ".");
        }
        if (line_number == 1) return;
        snplist.push_back(values[0]);
        double p_value = 1.0;
        if (values[6] != "." && values[6] != "NA" && values[6] != "NAN") p_value = atof(values[6].c_str());
        if (p_value < thresh && gws_snp_name_map.find(values[0]) == gws_snp_name_map.end()) {
            int size = gws_snp_name_map.size();
            gws_snp_name_map.insert(std::make_pair(values[0], size));
        }
    });
    return snplist;
}

double read_single_metafile(
    const string &filename,
    const map<string, int> &id_map,
    vector<string> &snp_a1,
    vector<string> &snp_a2,
    Eigen::VectorXd &snp_freq,
    Eigen::VectorXd &snp_b,
    Eigen::VectorXd &snp_se,
    Eigen::VectorXd &snp_pval,
    Eigen::VectorXd &snp_n,
    vector<bool> &snp_flag
) {
    vector<double> variances;
    for_each_line(filename, [&](int line_number, const string &line) {
        vector<string> values = fields(line);
        if (values.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + filename + "] should be in GCTA-COJO format, line " + std::to_string(line_number) + ".");
        }
        if (line_number == 1) return;
        auto iter = id_map.find(values[0]);
        if (iter == id_map.end()) return;
        int index = iter->second;
        bool missing = false;
        snp_a1[index] = uppercase(values[1]);
        snp_a2[index] = uppercase(values[2]);
        auto parse = [&](int column, Eigen::VectorXd &target) {
            const string &value = values[column];
            if (value != "." && value != "NA" && value != "NAN") target(index) = atof(value.c_str());
            else {
                target(index) = std::numeric_limits<double>::quiet_NaN();
                missing = true;
            }
        };
        parse(3, snp_freq);
        parse(4, snp_b);
        parse(5, snp_se);
        parse(6, snp_pval);
        parse(7, snp_n);
        snp_flag[index] = true;
        if (!missing) {
            double h = 2.0 * snp_freq(index) * (1.0 - snp_freq(index));
            variances.push_back(h * snp_b(index) * snp_b(index) + h * snp_n(index) * snp_se(index) * snp_se(index));
        }
    });
    return median(variances);
}

std::pair<double, double> prevalence(const vector<string> &values, const string &trait, const string &filename, int line_number) {
    double sample = std::numeric_limits<double>::quiet_NaN();
    double population = std::numeric_limits<double>::quiet_NaN();
    if (values.size() == 4) {
        string sample_text = uppercase(values[2]);
        string population_text = uppercase(values[3]);
        if (sample_text != "NA" && sample_text != "NAN" && sample_text != ".") {
            sample = atof(sample_text.c_str());
            if (sample < 0.0 || sample > 1.0) LOGGER.e(0, "invalid sample prevalence for trait [" + trait + "].");
        }
        if (population_text != "NA" && population_text != "NAN" && population_text != ".") {
            population = atof(population_text.c_str());
            if (population < 0.0 || population > 1.0) LOGGER.e(0, "invalid population prevalence for trait [" + trait + "].");
        }
    }
    (void)filename;
    (void)line_number;
    return std::make_pair(sample, population);
}

}

void SummaryData::init_meta_snp_map(vector<string> snplist, map<string, int> &snp_name_map, vector<string> &snp_name, vector<int> &remain_snp) {
    int size = 0;
    snp_name_map.clear();
    remain_snp.resize(snplist.size());
    for (size_t index = 0; index < snplist.size(); index++) {
        snp_name_map.insert(std::make_pair(snplist[index], index));
        if (size == snp_name_map.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[index] + ".");
        size = snp_name_map.size();
        remain_snp[index] = index;
    }
    snp_name = snplist;
}

void SummaryData::init_gwas_variable(vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, Eigen::MatrixXd &snp_freq, Eigen::MatrixXd &snp_b, Eigen::MatrixXd &snp_se, Eigen::MatrixXd &snp_pval, Eigen::MatrixXd &n, int trait_count, int snp_count) {
    snp_a1.clear();
    snp_a2.clear();
    snp_a1.resize(trait_count);
    snp_a2.resize(trait_count);
    for (int index = 0; index < trait_count; index++) {
        snp_a1[index].resize(snp_count);
        snp_a2[index].resize(snp_count);
    }
    snp_freq.resize(snp_count, trait_count);
    snp_b.resize(snp_count, trait_count);
    snp_se.resize(snp_count, trait_count);
    snp_pval.resize(snp_count, trait_count);
    n.resize(snp_count, trait_count);
}

void SummaryData::update_meta_snp_list(vector<string> &snplist, map<string, int> snp_id_map) {
    vector<string> original = snplist;
    snplist.clear();
    for (const string &name : original) {
        if (snp_id_map.find(name) != snp_id_map.end()) snplist.push_back(name);
    }
}

void SummaryData::update_meta_snp_map(vector<string> snplist, map<string, int> &snp_id_map, vector<string> &snp_id, vector<int> &snp_indx, bool indx_flag) {
    int size = 0;
    map<string, int> added;
    for (size_t index = 0, next = snp_id_map.size(); index < snplist.size(); index++) {
        added.insert(std::make_pair(snplist[index], index));
        if (size == added.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[index] + ".");
        size = added.size();
        if (snp_id_map.find(snplist[index]) != snp_id_map.end()) continue;
        snp_id_map.insert(std::make_pair(snplist[index], next));
        snp_id.push_back(snplist[index]);
        if (indx_flag) snp_indx.push_back(next);
        next++;
    }
}

void SummaryData::update_meta_snp(map<string, int> &snp_name_map, vector<string> &snp_name, vector<int> &snp_remain) {
    vector<string> old_names = snp_name;
    vector<int> old_remain = snp_remain;
    snp_name.resize(old_remain.size());
    snp_remain.resize(old_remain.size());
    for (size_t index = 0; index < old_remain.size(); index++) {
        snp_name[index] = old_names[old_remain[index]];
        snp_remain[index] = index;
    }
    snp_name_map.clear();
    for (size_t index = 0; index < snp_name.size(); index++) snp_name_map.insert(std::make_pair(snp_name[index], index));
}

void SummaryData::update_mtcojo_snp_rm(vector<string> adjusted_snps, map<string, int> &snp_id_map, vector<int> &remain_snp_indx) {
    for (const string &name : adjusted_snps) snp_id_map.erase(name);
    remain_snp_indx.clear();
    for (const auto &item : snp_id_map) remain_snp_indx.push_back(item.second);
    std::stable_sort(remain_snp_indx.begin(), remain_snp_indx.end());
}

vector<string> SummaryData::read_snp_metafile_txt(string metafile, map<string, int> &gws_snp_name_map, double thresh) {
    return read_snp_metafile(metafile, gws_snp_name_map, thresh);
}

vector<string> SummaryData::read_snp_metafile_gz(string metafile, map<string, int> &gws_snp_name_map, double thresh) {
    return read_snp_metafile(metafile, gws_snp_name_map, thresh);
}

double SummaryData::read_single_metafile_txt(string metafile, map<string, int> id_map, vector<string> &snp_a1, vector<string> &snp_a2, Eigen::VectorXd &snp_freq, Eigen::VectorXd &snp_b, Eigen::VectorXd &snp_se, Eigen::VectorXd &snp_pval, Eigen::VectorXd &snp_n, vector<bool> &snp_flag) {
    return read_single_metafile(metafile, id_map, snp_a1, snp_a2, snp_freq, snp_b, snp_se, snp_pval, snp_n, snp_flag);
}

double SummaryData::read_single_metafile_gz(string metafile, map<string, int> id_map, vector<string> &snp_a1, vector<string> &snp_a2, Eigen::VectorXd &snp_freq, Eigen::VectorXd &snp_b, Eigen::VectorXd &snp_se, Eigen::VectorXd &snp_pval, Eigen::VectorXd &snp_n, vector<bool> &snp_flag) {
    return read_single_metafile(metafile, id_map, snp_a1, snp_a2, snp_freq, snp_b, snp_se, snp_pval, snp_n, snp_flag);
}

void read_gsmr_file_list(string gsmr_file_list, vector<string> &pheno_name, vector<string> &pheno_file, vector<double> &popu_prev, vector<double> &smpl_prev) {
    std::ifstream input(gsmr_file_list);
    if (!input) LOGGER.e(0, "cannot open the file [" + gsmr_file_list + "] to read.");
    string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        line_number++;
        vector<string> values = fields(line);
        if (values.size() != 2 && values.size() != 4) {
            LOGGER.e(0, "the format of file [" + gsmr_file_list + "] is incorrect, line " + std::to_string(line_number) + ".");
        }
        pheno_name.push_back(values[0]);
        pheno_file.push_back(values[1]);
        auto prev = prevalence(values, pheno_name.back(), gsmr_file_list, line_number);
        smpl_prev.push_back(prev.first);
        popu_prev.push_back(prev.second);
    }
}

void read_metafile_list(string mtcojolist_file, string &target_pheno, string &target_pheno_file, vector<string> &covar_pheno, vector<string> &covar_pheno_file, vector<double> &popu_prev, vector<double> &smpl_prev) {
    std::ifstream input(mtcojolist_file);
    if (!input) LOGGER.e(0, "cannot open the file [" + mtcojolist_file + "] to read.");
    string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        line_number++;
        vector<string> values = fields(line);
        if (values.size() != 2 && values.size() != 4) {
            LOGGER.e(0, "the format of file [" + mtcojolist_file + "] is incorrect, line " + std::to_string(line_number) + ".");
        }
        auto prev = prevalence(values, values[0], mtcojolist_file, line_number);
        if (line_number == 1) {
            target_pheno = values[0];
            target_pheno_file = values[1];
        } else {
            covar_pheno.push_back(values[0]);
            covar_pheno_file.push_back(values[1]);
        }
        smpl_prev.push_back(prev.first);
        popu_prev.push_back(prev.second);
    }
}
