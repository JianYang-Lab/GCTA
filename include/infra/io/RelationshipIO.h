#ifndef GCTA2_RELATIONSHIPIO_H
#define GCTA2_RELATIONSHIPIO_H

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace RelationshipIO {

void write(
    const std::string &prefix,
    const std::vector<std::string> &family_ids,
    const std::vector<std::string> &individual_ids,
    const Eigen::MatrixXd &matrix,
    const Eigen::MatrixXd &denominators,
    bool binary
);

}

#endif
