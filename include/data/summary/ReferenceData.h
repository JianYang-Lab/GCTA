#ifndef GCTA2_REFERENCE_DATA_H
#define GCTA2_REFERENCE_DATA_H

#include "data/Geno.h"
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <string>
#include <vector>

class ReferenceData {
public:
    virtual ~ReferenceData() = default;

    std::map<std::string, int> _snp_name_map;
    std::vector<std::string> _snp_name;
    std::vector<int> _chr;
    std::vector<int> _bp;
    std::vector<std::string> _allele1;
    std::vector<std::string> _allele2;
    std::vector<std::string> _ref_A;
    std::vector<int> _include;
    std::vector<int> _keep;
    std::vector<double> _mu;

    void calcu_mu();
    bool make_XMat_subset(Eigen::MatrixXf &matrix, std::vector<int> &snp_indices, bool divide_by_std);
    void makex_eigenVector(int index, Eigen::VectorXd &values, bool resize = true, bool minus_2p = false);

protected:
    virtual void load() = 0;
    std::map<int, GenoBufItem> geno_items;
    bool loaded = false;
};

class V2ReferenceData : public ReferenceData {
public:
    V2ReferenceData();

private:
    std::unique_ptr<Pheno> pheno;
    std::unique_ptr<Marker> marker;
    std::unique_ptr<Geno> geno;
    void load() override;
};

#endif
