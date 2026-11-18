#ifndef GCTA2_SUMMARY_DATA_H
#define GCTA2_SUMMARY_DATA_H

#include <Eigen/Dense>

#include <map>
#include <string>
#include <vector>

class SummaryData {
public:
    int _expo_num = 0;
    int _outcome_num = 0;
    int _n_gsmr_rst_item = 5;
    int _gsmr_beta_version = 0;
    Eigen::MatrixXd _r_pheno_sample;
    std::vector<std::string> _gwas_trait_name;
    std::vector<std::vector<bool>> _snp_val_flag;

    std::string _target_pheno_name;
    std::vector<std::string> _meta_snp_name;
    std::vector<int> _meta_remain_snp;
    std::map<std::string, int> _meta_snp_name_map;
    std::vector<std::string> _covar_pheno_name;
    std::vector<std::string> _meta_snp_a1;
    std::vector<std::string> _meta_snp_a2;
    Eigen::MatrixXd _meta_snp_freq;
    Eigen::MatrixXd _meta_snp_b;
    Eigen::MatrixXd _meta_snp_se;
    Eigen::MatrixXd _meta_snp_pval;
    Eigen::MatrixXd _meta_snp_n_o;
    Eigen::VectorXd _meta_vp_trait;
    std::vector<double> _meta_popu_prev;
    std::vector<double> _meta_smpl_prev;

    void init_meta_snp_map(std::vector<std::string> snplist, std::map<std::string, int> &snp_name_map, std::vector<std::string> &snp_name, std::vector<int> &remain_snp);
    void init_gwas_variable(std::vector<std::vector<std::string>> &snp_a1, std::vector<std::vector<std::string>> &snp_a2, Eigen::MatrixXd &snp_freq, Eigen::MatrixXd &snp_b, Eigen::MatrixXd &snp_se, Eigen::MatrixXd &snp_pval, Eigen::MatrixXd &n, int trait_count, int snp_count);
    void update_meta_snp_list(std::vector<std::string> &snplist, std::map<std::string, int> snp_id_map);
    void update_meta_snp_map(std::vector<std::string> snplist, std::map<std::string, int> &snp_id_map, std::vector<std::string> &snp_id, std::vector<int> &snp_indx, bool indx_flag);
    void update_meta_snp(std::map<std::string, int> &snp_name_map, std::vector<std::string> &snp_name, std::vector<int> &snp_remain);
    void update_mtcojo_snp_rm(std::vector<std::string> adjusted_snps, std::map<std::string, int> &snp_id_map, std::vector<int> &remain_snp_indx);

    std::vector<std::string> read_snp_metafile_txt(std::string metafile, std::map<std::string, int> &gws_snp_name_map, double thresh);
    std::vector<std::string> read_snp_metafile_gz(std::string metafile, std::map<std::string, int> &gws_snp_name_map, double thresh);
    double read_single_metafile_txt(std::string metafile, std::map<std::string, int> id_map, std::vector<std::string> &snp_a1, std::vector<std::string> &snp_a2, Eigen::VectorXd &snp_freq, Eigen::VectorXd &snp_b, Eigen::VectorXd &snp_se, Eigen::VectorXd &snp_pval, Eigen::VectorXd &snp_n, std::vector<bool> &snp_flag);
    double read_single_metafile_gz(std::string metafile, std::map<std::string, int> id_map, std::vector<std::string> &snp_a1, std::vector<std::string> &snp_a2, Eigen::VectorXd &snp_freq, Eigen::VectorXd &snp_b, Eigen::VectorXd &snp_se, Eigen::VectorXd &snp_pval, Eigen::VectorXd &snp_n, std::vector<bool> &snp_flag);
};

void read_gsmr_file_list(std::string gsmr_file_list, std::vector<std::string> &pheno_name, std::vector<std::string> &pheno_file, std::vector<double> &popu_prev, std::vector<double> &smpl_prev);
void read_metafile_list(std::string mtcojolist_file, std::string &target_pheno, std::string &target_pheno_file, std::vector<std::string> &covar_pheno, std::vector<std::string> &covar_pheno_file, std::vector<double> &popu_prev, std::vector<double> &smpl_prev);

#endif
