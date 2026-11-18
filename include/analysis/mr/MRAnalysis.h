#ifndef GCTA2_MR_ANALYSIS_H
#define GCTA2_MR_ANALYSIS_H

#include "data/summary/ReferenceData.h"
#include "data/summary/SummaryData.h"
#include <Eigen/Dense>
#include <sstream>

using eigenMatrix = Eigen::MatrixXd;
using eigenVector = Eigen::VectorXd;

struct MROptions {
    std::string out;
    std::string exposure_list;
    std::string outcome_list;
    std::string mtcojo_list;
    std::string mtcojo_bxy_file;
    std::string ref_ld_dir;
    std::string w_ld_dir;
    double freq_thresh = 0.2;
    double gwas_thresh = 5e-8;
    double global_heidi_thresh = 0.01;
    double std_heidi_thresh = 0.01;
    double ld_fdr_thresh = 0.05;
    double clump_wind_size = 10000;
    double clump_r2_thresh = 0.05;
    int nsnp_gsmr = 10;
    int gsmr_direction = 0;
    int gsmr_so_alg = -9;
    int gsmr_beta_version = 0;
    bool effect_plot = false;
};

MROptions consume_common_mr_options(std::map<std::string, std::vector<std::string>> &options_in);

namespace MR {
void update_id_map_kp(const std::vector<std::string> &id_list, std::map<std::string, int> &id_map, std::vector<int> &keep);
void update_id_map_rm(const std::vector<std::string> &id_list, std::map<std::string, int> &id_map, std::vector<int> &keep);
}

class MRAnalysis : public SummaryData, public V2ReferenceData {
public:
    explicit MRAnalysis(std::string out) : _out(std::move(out)) {}

    std::string _out;

    void read_gsmrfile(std::string expo_file_list, std::string outcome_file_list, double gwas_thresh, int nsnp_gsmr, int gsmr_so_alg);
    void gsmr(int gsmr_alg_flag, std::string ref_ld_dirt, std::string w_ld_dirt, double freq_thresh, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, bool o_snp_instru_flag, int gsmr_so_alg, int gsmr_beta_version);
    std::vector<std::vector<double>> forward_gsmr(std::stringstream &ss, std::map<std::string, int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, std::stringstream &ss_pleio);
    std::vector<std::vector<double>> reverse_gsmr(std::stringstream &ss, std::map<std::string, int> &snp_instru_map, double gwas_thresh, double clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, std::stringstream &ss_pleio);
    eigenMatrix rho_sample_overlap(std::vector<std::vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, std::vector<std::string> snp_name, std::vector<int> snp_remain, std::string ref_ld_dirt, std::string w_ld_dirt, std::vector<std::string> trait_name, int gsmr_so_alg);
    eigenMatrix sample_overlap_ldsc(std::vector<std::vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int nexpo, int noutcome, std::vector<std::string> snp_name, std::vector<int> snp_remain, std::string ref_ld_dirt, std::string w_ld_dirt, std::vector<std::string> trait_name);
    eigenMatrix sample_overlap_rb(std::vector<std::vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, int nexpo, int noutcome, std::vector<std::string> snp_name, std::vector<int> snp_remain, std::vector<std::string> trait_name);

    int read_mtcojofile(std::string mtcojolist_file, double gwas_thresh, int nsnp_gsmr);
    void mtcojo(std::string mtcojo_bxy_file, std::string ref_ld_dirt, std::string w_ld_dirt, double freq_thresh, double gwas_thresh, int clump_wind_size, double clump_r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, int gsmr_beta_version);

    std::vector<std::string> remove_bad_snps(std::vector<std::string> snp_name, std::vector<int> snp_remain, std::vector<std::vector<bool>> snp_flag, std::vector<std::vector<std::string>> &snp_a1, std::vector<std::vector<std::string>> &snp_a2, eigenMatrix &snp_freq, eigenMatrix &snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, std::map<std::string, int> plink_snp_name_map, std::vector<std::string> snp_ref_a1, std::vector<std::string> snp_ref_a2, int ntarget, int ncovar, std::string outfile_name);
    std::vector<std::string> remove_freq_diff_snps(std::vector<std::string> meta_snp_name, std::vector<int> meta_snp_remain, std::map<std::string, int> snp_name_map, std::vector<double> ref_freq, eigenMatrix meta_freq, std::vector<std::vector<bool>> snp_flag, int ntrait, double freq_thresh, std::string outfile_name);
    std::vector<std::string> remove_mono_snps(std::map<std::string, int> snp_name_map, std::vector<double> ref_snpfreq, std::string outfile_name);
    std::vector<std::string> filter_meta_snp_pval(std::vector<std::string> snp_name, std::vector<int> remain_snp_indx, eigenMatrix snp_pval, int start_indx, int end_indx, std::vector<std::vector<bool>> snp_flag, double pval_thresh);

    std::vector<std::string> clumping_meta(eigenVector snp_chival, std::vector<bool> snp_flag, double pval_thresh, int wind_size, double r2_thresh);
    std::vector<double> gsmr_meta(std::vector<std::string> &snp_instru, eigenVector bzx, eigenVector bzx_se, eigenVector bzx_pval, eigenVector bzy, eigenVector bzy_se, eigenVector bzy_pval, double rho_pheno, std::vector<bool> snp_flag, double gwas_thresh, int wind_size, double r2_thresh, double std_heidi_thresh, double global_heidi_thresh, double ld_fdr_thresh, int nsnp_gsmr, std::string &pleio_snps, std::string &err_msg);

    std::vector<std::string> read_snp_ldsc(std::map<std::string, int> ldsc_snp_name_map, std::vector<std::string> snp_name, std::vector<int> snp_remain, int &ttl_mk_num, std::string ref_ld_dirt, std::string w_ld_dirt, std::vector<double> &ref_ld_vec, std::vector<double> &w_ld_vec);
    void reorder_snp_effect(std::vector<int> snp_remain, eigenMatrix &bhat_z, eigenMatrix &bhat_n, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, std::vector<std::vector<bool>> &snp_flag, std::vector<std::vector<bool>> snp_val_flag, std::vector<int> &nsnp_cm_trait, std::vector<std::string> cm_ld_snps, std::map<std::string, int> ldsc_snp_name_map, eigenVector &ref_ld, eigenVector &w_ld, std::vector<double> ref_ld_vec, std::vector<double> w_ld_vec, int ntrait);
    eigenMatrix ldsc_snp_h2(eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, std::vector<std::vector<bool>> snp_flag, std::vector<int> nsnp_cm_trait, int n_cm_ld_snps, int ttl_mk_num, std::vector<std::string> trait_name, int ntrait);
    eigenMatrix ldsc_snp_rg(eigenMatrix ldsc_var_h2, eigenMatrix bhat_z, eigenMatrix bhat_n, eigenVector ref_ld, eigenVector w_ld, std::vector<std::vector<bool>> snp_flag, std::vector<int> trait_indx1, std::vector<int> trait_indx2, int n_cm_ld_snps, int ttl_mk_num, std::vector<std::string> trait_name);
    bool mtcojo_ldsc(std::vector<std::vector<bool>> snp_val_flag, eigenMatrix snp_b, eigenMatrix snp_se, eigenMatrix snp_n, int ntrait, std::vector<std::string> snp_name, std::vector<int> snp_remain, std::string ref_ld_dirt, std::string w_ld_dirt, std::vector<std::string> trait_name, eigenMatrix &ldsc_intercept, eigenMatrix &ldsc_slope);
};

#endif
