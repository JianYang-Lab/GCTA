#include "analysis/mr/MRAnalysis.h"

#include "infrastructure/Logger.h"
#include <algorithm>
#include <cstdlib>

using std::map;
using std::string;
using std::vector;

namespace {

vector<string> take(map<string, vector<string>> &options, const string &flag) {
    vector<string> values;
    auto iter = options.find(flag);
    if (iter == options.end()) return values;
    values = iter->second;
    options.erase(iter);
    return values;
}

double number(const vector<string> &values, double current) {
    return values.empty() ? current : atof(values[0].c_str());
}

}

MROptions consume_common_mr_options(map<string, vector<string>> &options_in) {
    MROptions options;
    options.out = options_in["out"][0];
    auto values = take(options_in, "--ref-ld-chr");
    if (!values.empty()) options.ref_ld_dir = values[0];
    values = take(options_in, "--w-ld-chr");
    if (!values.empty()) options.w_ld_dir = values[0];
    if (options.ref_ld_dir.empty() != options.w_ld_dir.empty()) {
        if (options.ref_ld_dir.empty()) LOGGER.e(0, "--w-ld-chr, please specify the directory of the LD scores for the regression weights.");
        LOGGER.e(0, "--ref-ld-chr, please specify the directory of the LD score files.");
    }
    if (!options.ref_ld_dir.empty() && options.ref_ld_dir.back() != '/') options.ref_ld_dir += '/';
    if (!options.w_ld_dir.empty() && options.w_ld_dir.back() != '/') options.w_ld_dir += '/';

    options.freq_thresh = number(take(options_in, "--diff-freq"), options.freq_thresh);
    if (options.freq_thresh < 0.0 || options.freq_thresh > 1.0) LOGGER.e(0, "--diff-freq, Invalid threshold used to check allele frequency difference.");
    options.gwas_thresh = number(take(options_in, "--gwas-thresh"), options.gwas_thresh);
    if (options.gwas_thresh < 0.0 || options.gwas_thresh > 1.0) LOGGER.e(0, "--gwas-thresh, Invalid GWAS p-value threshold.");
    values = take(options_in, "--heidi-thresh");
    if (!values.empty()) {
        options.std_heidi_thresh = atof(values[0].c_str());
        if (values.size() > 1) options.global_heidi_thresh = atof(values[1].c_str());
        else if (options.std_heidi_thresh == 0.0) options.global_heidi_thresh = 0.0;
        if (values.size() > 2) LOGGER.e(0, "--heidi-thresh, please specify p-value threshold(s) for the HEIDI-outlier analysis.");
    }
    if (options.std_heidi_thresh < 0.0 || options.std_heidi_thresh > 1.0) LOGGER.e(0, "--heidi-thresh, Invalid p-value threshold for single-SNP-based HEIDI-outlier test.");
    if (options.global_heidi_thresh < 0.0 || options.global_heidi_thresh > 1.0) LOGGER.e(0, "--heidi-thresh, Invalid p-value threshold for multi-SNP-based HEIDI-outlier test.");
    options.ld_fdr_thresh = number(take(options_in, "--gsmr-ld-fdr"), options.ld_fdr_thresh);
    if (options.ld_fdr_thresh < 0.0 || options.ld_fdr_thresh > 1.0) LOGGER.e(0, "--gsmr-ld-fdr, Invalid FDR threshold for LD correlation matrix.");
    options.clump_wind_size = number(take(options_in, "--clump-kb"), options.clump_wind_size);
    if (options.clump_wind_size < 0.0 || options.clump_wind_size > 1e6) LOGGER.e(0, "--clump-kb, Invalid window size for the clumping analysis.");
    options.clump_r2_thresh = number(take(options_in, "--clump-r2"), options.clump_r2_thresh);
    if (options.clump_r2_thresh < 0.0 || options.clump_r2_thresh > 1.0) LOGGER.e(0, "--clump-r2, Invalid LD r2 threshold for the clumping analysis.");

    values = take(options_in, "--gsmr-snp-min");
    if (values.empty()) values = take(options_in, "--gsmr-snp");
    if (!values.empty()) options.nsnp_gsmr = atoi(values[0].c_str());
    if (options.nsnp_gsmr < 0 || options.nsnp_gsmr > 1e6) LOGGER.e(0, "--gsmr-snp-min, Invalid SNP number threshold for the GSMR analysis.");
    values = take(options_in, "--gsmr-direction");
    if (!values.empty()) options.gsmr_direction = atoi(values[0].c_str());
    if (options.gsmr_direction < 0 || options.gsmr_direction > 2) LOGGER.e(0, "--gsmr-direction should be 0 (forward-GSMR), 1 (reverse-GSMR) or 2 (bi-GSMR).");
    values = take(options_in, "--gsmr-so");
    if (!values.empty()) options.gsmr_so_alg = atoi(values[0].c_str());
    if (options.gsmr_so_alg != -9 && options.gsmr_so_alg != 0 && options.gsmr_so_alg != 1) {
        LOGGER.e(0, "--gsmr-so should be 0 (LD score regression) or 1 (correlation of SNP effects).");
    }
    if (options_in.find("--gsmr2-beta") != options_in.end()) {
        options.gsmr_beta_version = 1;
        options_in.erase("--gsmr2-beta");
    }
    if (options_in.find("--effect-plot") != options_in.end()) {
        options.effect_plot = true;
        options_in.erase("--effect-plot");
    }
    if (options_in.find("--gsmr-alg") != options_in.end()) LOGGER.e(0, "--gsmr-alg has been superseded by --gsmr-direction.");
    if (options_in.find("--heidi-snp") != options_in.end()) LOGGER.e(0, "--heidi-snp is discontinued. Please use --gsmr-snp-min to specify minimum number of SNP instruments for the HEIDI-outlier analysis.");
    if (options_in.find("--clump-p1") != options_in.end()) LOGGER.e(0, "--clump-p1 is discontinued. Please use --gwas-thresh to specify a p-value threshold for index SNPs.");
    return options;
}

namespace MR {

void update_id_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep) {
    map<string, int> remaining(id_map);
    for (const string &id : id_list) remaining.erase(id);
    for (const auto &item : remaining) id_map.erase(item.first);
    keep.clear();
    for (const auto &item : id_map) keep.push_back(item.second);
    std::stable_sort(keep.begin(), keep.end());
}

void update_id_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep) {
    for (const string &id : id_list) id_map.erase(id);
    keep.clear();
    for (const auto &item : id_map) keep.push_back(item.second);
    std::stable_sort(keep.begin(), keep.end());
}

}
