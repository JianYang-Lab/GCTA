#ifndef GCTA2_ERM_H
#define GCTA2_ERM_H

#include "data/expression/ExpressionData.h"
#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

class ERM {
public:
    ERM(const ExpressionData &data, int algorithm);

    static int registerOption(std::map<std::string, std::vector<std::string>> &options_in);
    static void processMain();
    void write(const std::string &prefix, bool binary) const;

private:
    const ExpressionData &data;
    int algorithm;
    Eigen::MatrixXd transformed;
    Eigen::MatrixXd matrix;
    Eigen::MatrixXd denominators;

    void transform();
    void calculate();
};

#endif
