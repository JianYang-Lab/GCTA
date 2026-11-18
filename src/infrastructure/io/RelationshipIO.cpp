#include "infrastructure/io/RelationshipIO.h"

#include "infrastructure/Logger.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <zlib.h>

using std::string;
using std::vector;

namespace RelationshipIO {

static void write_ids(
    const string &prefix,
    const vector<string> &family_ids,
    const vector<string> &individual_ids,
    const string &matrix_file
) {
    string filename = prefix + ".grm.id";
    std::ofstream output(filename);
    if (!output) LOGGER.e(0, "cannot open file [" + filename + "] to write.");
    for (size_t index = 0; index < family_ids.size(); index++) {
        output << family_ids[index] << "\t" << individual_ids[index] << std::endl;
    }
    LOGGER << "IDs for the GRM file [" + matrix_file + "] have been saved in the file [" + filename + "]." << std::endl;
}

static string write_binary(
    const string &prefix,
    const Eigen::MatrixXd &matrix,
    const Eigen::MatrixXd &denominators
) {
    string matrix_file = prefix + ".grm.bin";
    std::fstream matrix_output(matrix_file, std::ios::out | std::ios::binary);
    if (!matrix_output) LOGGER.e(0, "cannot open file [" + matrix_file + "] to write.");
    string count_file = prefix + ".grm.N.bin";
    std::fstream count_output(count_file, std::ios::out | std::ios::binary);
    if (!count_output) LOGGER.e(0, "cannot open file [" + count_file + "] to write.");

    int sample_count = matrix.rows();
    for (int row = 0; row < sample_count; row++) {
        for (int column = 0; column <= row; column++) {
            float value = static_cast<float>(matrix(row, column));
            matrix_output.write(reinterpret_cast<char *>(&value), sizeof(value));
            float denominator = static_cast<float>(denominators(row, column));
            count_output.write(reinterpret_cast<char *>(&denominator), sizeof(denominator));
        }
    }
    if (!matrix_output.good() || !count_output.good()) {
        LOGGER.e(0, "failed to write the relationship matrix [" + prefix + "].");
    }
    matrix_output.close();
    count_output.close();
    LOGGER << "GRM of " << sample_count << " individuals has been saved in the file [" + matrix_file + "] (in binary format)." << std::endl;
    LOGGER << "Number of SNPs to calculate the genetic relationship between each pair of individuals has been saved in the file [" + count_file + "] (in binary format)." << std::endl;
    return matrix_file;
}

static string write_gzip(
    const string &prefix,
    const Eigen::MatrixXd &matrix,
    const Eigen::MatrixXd &denominators
) {
    string matrix_file = prefix + ".grm.gz";
    gzFile output = gzopen(matrix_file.c_str(), "wb");
    if (!output) LOGGER.e(0, "cannot open file [" + matrix_file + "] to write.");
    LOGGER << "Saving the genetic relationship matrix to the file [" + matrix_file + "] (in compressed text format)." << std::endl;
    for (int row = 0; row < matrix.rows(); row++) {
        for (int column = 0; column <= row; column++) {
            std::ostringstream line;
            line << std::scientific << std::setprecision(6)
                 << row + 1 << '\t' << column + 1 << '\t'
                 << denominators(row, column) << '\t' << matrix(row, column) << std::endl;
            string data = line.str();
            if (gzwrite(output, data.data(), static_cast<unsigned int>(data.size())) != static_cast<int>(data.size())) {
                gzclose(output);
                LOGGER.e(0, "failed to write [" + matrix_file + "].");
            }
        }
    }
    if (gzclose(output) != Z_OK) LOGGER.e(0, "failed to close [" + matrix_file + "].");
    LOGGER << "The genetic relationship matrix has been saved in the file [" + matrix_file + "] (in compressed text format)." << std::endl;
    return matrix_file;
}

void write(
    const string &prefix,
    const vector<string> &family_ids,
    const vector<string> &individual_ids,
    const Eigen::MatrixXd &matrix,
    const Eigen::MatrixXd &denominators,
    bool binary
) {
    if (family_ids.size() != individual_ids.size() || family_ids.size() != static_cast<size_t>(matrix.rows())) {
        LOGGER.e(0, "inconsistent relationship matrix and sample IDs.");
    }
    string matrix_file = binary ? write_binary(prefix, matrix, denominators) : write_gzip(prefix, matrix, denominators);
    write_ids(prefix, family_ids, individual_ids, matrix_file);
}

}
