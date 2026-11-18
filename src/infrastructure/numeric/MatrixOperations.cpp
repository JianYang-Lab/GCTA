#include "infrastructure/numeric/MatrixOperations.h"

#include "infrastructure/cpu_f77blas.h"
#include "infrastructure/Logger.h"
#include <vector>

namespace MatrixOperations {

bool inverse_llt(Eigen::MatrixXd &matrix, double &logdet) {
    int size = matrix.cols();
    std::vector<double> values(size * size);
    for (int row = 0; row < size; row++) {
        for (int column = 0; column < size; column++) values[row * size + column] = matrix(row, column);
    }

    int info = 0;
    char uplo = 'L';
#if GCTA_CPU_x86
    dpotrf(&uplo, &size, values.data(), &size, &info);
#else
    dpotrf_(&uplo, &size, values.data(), &size, &info);
#endif
    if (info < 0) LOGGER.e(0, "Cholesky decomposition failed. Invalid values found in the matrix.\n");
    if (info > 0) return false;

    logdet = 0.0;
    for (int index = 0; index < size; index++) {
        double diagonal = values[index * size + index];
        logdet += log(diagonal * diagonal);
    }
#if GCTA_CPU_x86
    dpotri(&uplo, &size, values.data(), &size, &info);
#else
    dpotri_(&uplo, &size, values.data(), &size, &info);
#endif
    if (info < 0) LOGGER.e(0, "invalid values found in the variance-covariance (V) matrix.\n");
    if (info > 0) return false;

    for (int column = 0; column < size; column++) {
        for (int row = 0; row <= column; row++) {
            matrix(row, column) = matrix(column, row) = values[row * size + column];
        }
    }
    return true;
}

}
