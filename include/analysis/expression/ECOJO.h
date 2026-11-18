#ifndef GCTA2_ECOJO_H
#define GCTA2_ECOJO_H

#include "data/expression/ExpressionData.h"

#include <Eigen/Dense>

#include <map>
#include <string>
#include <vector>

class ECOJO {
public:
    explicit ECOJO(std::string output);

    static int registerOption(std::map<std::string, std::vector<std::string>> &options_in);
    static void processMain();

    void load_expression(const ExpressionData &data);
    void load_correlation(const std::string &filename);
    void read_summary(const std::string &filename);
    void select(double p_cutoff, double collinear);
    void blup(double lambda);

private:
    std::string output;
    const ExpressionData *expression_data = nullptr;
    std::vector<std::string> probe_names;
    std::map<std::string, int> probe_map;
    std::vector<int> include;
    size_t sample_count = 0;
    Eigen::MatrixXd whole_correlation;
    Eigen::MatrixXd correlation;
    Eigen::MatrixXd cross_correlation;
    Eigen::VectorXd z_values;
    Eigen::VectorXd b_values;
    Eigen::VectorXd se_values;
    Eigen::VectorXd n_values;
    Eigen::VectorXd p_values;
    double p_cutoff = 5e-6;
    double collinear = 0.9;

    void calculate_correlation(const ExpressionData &data);
    void joint(const std::vector<int> &selected, Eigen::VectorXd &b_joint, Eigen::VectorXd &se_joint, Eigen::VectorXd &p_joint);
    void conditional(const std::vector<int> &selected, const std::vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond);
    void stepwise(std::vector<int> &selected, std::vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond);
    bool select_entry(std::vector<int> &selected, std::vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond);
    void select_stay(std::vector<int> &selected, Eigen::VectorXd &b_joint, Eigen::VectorXd &se_joint, Eigen::VectorXd &p_joint);
    bool init_correlation(const std::vector<int> &selected);
    void init_cross_correlation(const std::vector<int> &selected, const std::vector<int> &remain);
    bool insert_correlation(const std::vector<int> &selected, int insert_index);
    void erase_correlation(const std::vector<int> &selected);
    void invert_correlation();
    void write_selection(const std::vector<int> &selected, const Eigen::VectorXd &b_joint, const Eigen::VectorXd &se_joint, const Eigen::VectorXd &p_joint);
};

#endif
