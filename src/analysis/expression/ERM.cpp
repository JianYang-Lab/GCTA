#include "analysis/expression/ERM.h"

#include "infra/io/RelationshipIO.h"
#include "infra/Logger.h"

using std::map;
using std::string;
using std::vector;

namespace {
map<string, string> options;
bool output_binary = true;
}

ERM::ERM(const ExpressionData &data, int algorithm) : data(data), algorithm(algorithm) {
    transform();
    calculate();
}

void ERM::transform() {
    transformed = data.values();
    int sample_count = data.sample_count();
    int probe_count = data.probe_count();
    const vector<vector<bool>> &valid = data.valid();

    if (algorithm < 2) {
        Eigen::VectorXd means = Eigen::VectorXd::Zero(probe_count);
        Eigen::VectorXd nonmissing = Eigen::VectorXd::Zero(probe_count);
        for (int probe = 0; probe < probe_count; probe++) {
            for (int sample = 0; sample < sample_count; sample++) {
                if (valid[sample][probe]) {
                    means(probe) += transformed(sample, probe);
                    nonmissing(probe) += 1.0;
                }
            }
            means(probe) /= nonmissing(probe);
        }
        for (int sample = 0; sample < sample_count; sample++) {
            for (int probe = 0; probe < probe_count; probe++) {
                transformed(sample, probe) = valid[sample][probe]
                    ? transformed(sample, probe) - means(probe)
                    : 0.0;
            }
        }
        if (algorithm == 0) {
            Eigen::VectorXd scales(probe_count);
            for (int probe = 0; probe < probe_count; probe++) {
                scales(probe) = sqrt(transformed.col(probe).dot(transformed.col(probe)) / (nonmissing(probe) - 1.0));
            }
            for (int sample = 0; sample < sample_count; sample++) {
                for (int probe = 0; probe < probe_count; probe++) {
                    transformed(sample, probe) /= scales(probe);
                }
            }
        }
    } else {
        Eigen::VectorXd means = Eigen::VectorXd::Zero(sample_count);
        Eigen::VectorXd nonmissing = Eigen::VectorXd::Zero(sample_count);
        for (int sample = 0; sample < sample_count; sample++) {
            for (int probe = 0; probe < probe_count; probe++) {
                if (valid[sample][probe]) {
                    means(sample) += transformed(sample, probe);
                    nonmissing(sample) += 1.0;
                }
            }
            means(sample) /= nonmissing(sample);
        }
        for (int sample = 0; sample < sample_count; sample++) {
            for (int probe = 0; probe < probe_count; probe++) {
                transformed(sample, probe) = valid[sample][probe]
                    ? transformed(sample, probe) - means(sample)
                    : 0.0;
            }
        }
    }
}

void ERM::calculate() {
    LOGGER << "Recoding gene expression / methylation data ..." << std::endl;
    int sample_count = data.sample_count();
    int probe_count = data.probe_count();
    const vector<vector<bool>> &valid = data.valid();
    vector<vector<int>> missing(sample_count);
    for (int sample = 0; sample < sample_count; sample++) {
        for (int probe = 0; probe < probe_count; probe++) {
            if (!valid[sample][probe]) missing[sample].push_back(probe);
        }
    }

    denominators.setZero(sample_count, sample_count);
    if (algorithm == 0) {
        #pragma omp parallel for
        for (int row = 0; row < sample_count; row++) {
            for (int column = 0; column <= row; column++) {
                int shared_missing = 0;
                for (int probe : missing[column]) shared_missing += static_cast<int>(valid[row][probe]);
                denominators(row, column) = probe_count - missing[row].size() - shared_missing;
            }
        }
    } else if (algorithm == 1) {
        Eigen::VectorXd nonmissing = Eigen::VectorXd::Zero(probe_count);
        for (int probe = 0; probe < probe_count; probe++) {
            for (int sample = 0; sample < sample_count; sample++) {
                nonmissing(probe) += static_cast<int>(valid[sample][probe]);
            }
        }
        Eigen::VectorXd variances(probe_count);
        for (int probe = 0; probe < probe_count; probe++) {
            variances(probe) = nonmissing(probe) > 1.0
                ? transformed.col(probe).dot(transformed.col(probe)) / (nonmissing(probe) - 1.0)
                : 0.0;
        }
        double sum_variances = variances.sum();
        #pragma omp parallel for
        for (int row = 0; row < sample_count; row++) {
            double row_missing = 0.0;
            for (int probe : missing[row]) row_missing += variances(probe);
            for (int column = 0; column <= row; column++) {
                double column_missing = 0.0;
                for (int probe : missing[column]) {
                    if (valid[row][probe]) column_missing += variances(probe);
                }
                denominators(row, column) = sum_variances - row_missing - column_missing;
            }
        }
    } else {
        Eigen::VectorXd sum_squares(sample_count);
        for (int sample = 0; sample < sample_count; sample++) {
            sum_squares(sample) = transformed.row(sample).dot(transformed.row(sample));
        }
        #pragma omp parallel for
        for (int row = 0; row < sample_count; row++) {
            for (int column = 0; column <= row; column++) {
                double row_ssx = sum_squares(row);
                double column_ssx = sum_squares(column);
                for (int probe : missing[column]) {
                    if (valid[row][probe]) row_ssx -= transformed(row, probe) * transformed(row, probe);
                }
                for (int probe : missing[row]) {
                    if (valid[column][probe]) column_ssx -= transformed(column, probe) * transformed(column, probe);
                }
                denominators(row, column) = sqrt(row_ssx * column_ssx);
            }
        }
    }

    LOGGER << "\nCalculating expression relationship matrix (ERM) ... " << std::endl;
    matrix = transformed * transformed.transpose();
    for (int row = 0; row < sample_count; row++) {
        for (int column = 0; column <= row; column++) {
            matrix(row, column) = denominators(row, column) > 0.0
                ? matrix(row, column) / denominators(row, column)
                : 0.0;
        }
    }
    matrix = matrix.array() / matrix.diagonal().mean();
}

void ERM::write(const string &prefix, bool binary) const {
    RelationshipIO::write(prefix, data.family_ids(), data.individual_ids(), matrix, denominators, binary);
}

int ERM::registerOption(map<string, vector<string>> &options_in) {
    int result = 0;
    options["out"] = options_in["out"][0];
    const string expression_flag = "--efile";
    const string erm_flag = "--make-erm";
    if (options_in.find(erm_flag) != options_in.end()) {
        options["algorithm"] = "1";
        output_binary = true;
        result++;
        options_in.erase(erm_flag);
    }
    const string gz_flag = "--make-erm-gz";
    if (options_in.find(gz_flag) != options_in.end()) {
        options["algorithm"] = "1";
        output_binary = false;
        result++;
        options_in.erase(gz_flag);
    }
    const string algorithm_flag = "--make-erm-alg";
    if (options_in.find(algorithm_flag) != options_in.end()) {
        if (options_in[algorithm_flag].size() == 1) options["algorithm"] = options_in[algorithm_flag][0];
        result++;
        options_in.erase(algorithm_flag);
    }
    if (result > 0 && options_in.find(expression_flag) != options_in.end()) {
        options["efile"] = options_in[expression_flag][0];
        options_in.erase(expression_flag);
    }
    return result;
}

void ERM::processMain() {
    if (options.find("efile") == options.end()) LOGGER.e(0, "please specify the expression data with --efile.");
    int algorithm = 0;
    try {
        if (options.find("algorithm") != options.end()) algorithm = std::stoi(options["algorithm"]);
    } catch (std::invalid_argument &) {
        LOGGER.e(0, "--make-erm-alg should be 1, 2 or 3.");
    }
    if (algorithm < 1 || algorithm > 3) LOGGER.e(0, "\n --make-erm-alg should be 1, 2 or 3.\n");
    algorithm--;
    ExpressionData data(options["efile"]);
    ERM erm(data, algorithm);
    erm.write(options["out"] + ".E", output_binary);
}
