#ifndef GCTA2_MATRIXOPERATIONS_H
#define GCTA2_MATRIXOPERATIONS_H

#include <Eigen/Dense>

namespace MatrixOperations {

bool inverse_llt(Eigen::MatrixXd &matrix, double &logdet);

}

#endif
