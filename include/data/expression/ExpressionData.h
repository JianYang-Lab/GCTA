#ifndef GCTA2_EXPRESSIONDATA_H
#define GCTA2_EXPRESSIONDATA_H

#include <Eigen/Dense>

#include <string>
#include <vector>

class ExpressionData {
public:
    explicit ExpressionData(const std::string &filename);

    size_t sample_count() const { return fid.size(); }
    size_t probe_count() const { return probe_names.size(); }
    const std::vector<std::string> &family_ids() const { return fid; }
    const std::vector<std::string> &individual_ids() const { return pid; }
    const std::vector<std::string> &probes() const { return probe_names; }
    const Eigen::MatrixXd &values() const { return probe_data; }
    const std::vector<std::vector<bool>> &valid() const { return valid_data; }

private:
    std::vector<std::string> fid;
    std::vector<std::string> pid;
    std::vector<std::string> probe_names;
    Eigen::MatrixXd probe_data;
    std::vector<std::vector<bool>> valid_data;
};

#endif
