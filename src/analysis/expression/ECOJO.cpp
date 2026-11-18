#include "analysis/expression/ECOJO.h"

#include "analysis/StatLib.h"

#include "infra/Logger.h"
#include "infra/numeric/MatrixOperations.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <utility>

using std::map;
using std::string;
using std::vector;

namespace {
map<string, string> options;
string action;
}

ECOJO::ECOJO(string output) : output(std::move(output)) {}

void ECOJO::load_expression(const ExpressionData &data) {
    expression_data = &data;
    sample_count = data.sample_count();
    probe_names = data.probes();
    for (size_t index = 0; index < probe_names.size(); index++) {
        probe_map.insert(std::make_pair(probe_names[index], index));
        include.push_back(index);
    }
}

void ECOJO::load_correlation(const string &filename) {
    std::ifstream input(filename);
    if (!input.is_open()) LOGGER.e(0, "cannot open the file [" + filename + "] to read.");
    LOGGER << "Reading correlation matrix of gene expression from [" + filename + "] ..." << std::endl;
    string line;
    std::getline(input, line);
    std::stringstream header(line);
    string name;
    while (header >> name) {
        probe_map.insert(std::make_pair(name, static_cast<int>(probe_names.size())));
        if (probe_map.size() != probe_names.size() + 1) LOGGER.e(0, "duplicated probe names found: \"" + name + "\".");
        probe_names.push_back(name);
        include.push_back(include.size());
    }
    LOGGER << probe_names.size() << " probes found in the file. \nReading correlation matrix ..." << std::endl;
    whole_correlation.resize(probe_names.size(), probe_names.size());
    for (size_t row = 0; row < probe_names.size(); row++) {
        for (size_t column = 0; column < probe_names.size(); column++) {
            if (!(input >> whole_correlation(row, column))) LOGGER.e(0, "incorrect format of [" + filename + "].");
        }
    }
    LOGGER << "Correlation matrix for " << probe_names.size() << " probes have been included from the file [" + filename + "]." << std::endl;
}

void ECOJO::read_summary(const string &filename) {
    LOGGER << "\nReading expression-trait association summary-level statistics from [" + filename + "] ..." << std::endl;
    std::ifstream input(filename);
    if (!input) LOGGER.e(0, "cannot open the file [" + filename + "] to read.");
    string line;
    std::getline(input, line);
    vector<string> header;
    std::stringstream header_stream(line);
    string item;
    while (header_stream >> item) header.push_back(item);
    if (header.size() < 3) LOGGER.e(0, "there needs to be at least 3 columns in the file [" + filename + "].");

    map<string, std::pair<double, double>> summary;
    int matched_rows = 0;
    int line_number = 1;
    while (std::getline(input, line)) {
        line_number++;
        std::stringstream stream(line);
        string probe;
        double z_value, n_value;
        if (!(stream >> probe)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
        if (probe_map.find(probe) == probe_map.end()) continue;
        if (!(stream >> z_value)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
        if (!(stream >> n_value)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
        summary.insert(std::make_pair(probe, std::make_pair(z_value, n_value)));
        matched_rows++;
    }
    if (matched_rows == 0) LOGGER.e(0, "no probe remains in the analysis.");
    LOGGER << "GWAS summary statistics of " << matched_rows << " probs read from [" + filename + "]." << std::endl;
    LOGGER << "Matching the summary data to the genotype data ..." << std::endl;

    include.clear();
    for (size_t index = 0; index < probe_names.size(); index++) {
        if (summary.find(probe_names[index]) != summary.end()) include.push_back(index);
    }
    int count = include.size();
    z_values.resize(count);
    b_values.resize(count);
    se_values.resize(count);
    n_values.resize(count);
    p_values.resize(count);
    for (int index = 0; index < count; index++) {
        const auto &values = summary[probe_names[include[index]]];
        z_values(index) = values.first;
        n_values(index) = values.second;
        p_values(index) = StatLib::pchisq(z_values(index) * z_values(index), 1);
        double scale = sqrt(z_values(index) * z_values(index) + n_values(index));
        if (scale < 1e-30) {
            b_values(index) = 0.0;
            se_values(index) = -1.0;
        } else {
            b_values(index) = z_values(index) / scale;
            se_values(index) = 1.0 / scale;
        }
    }
    if (expression_data) calculate_correlation(*expression_data);
}

void ECOJO::calculate_correlation(const ExpressionData &data) {
    int probe_count = include.size();
    int individual_count = sample_count;
    Eigen::MatrixXd centered(probe_count, individual_count);
    const vector<vector<bool>> &valid = data.valid();
    for (int probe = 0; probe < probe_count; probe++) {
        double mean = 0.0;
        double nonmissing = 0.0;
        int source_probe = include[probe];
        for (int sample = 0; sample < individual_count; sample++) {
            if (valid[sample][source_probe]) {
                mean += data.values()(sample, source_probe);
                nonmissing += 1.0;
            }
        }
        mean /= nonmissing;
        for (int sample = 0; sample < individual_count; sample++) {
            centered(probe, sample) = valid[sample][source_probe]
                ? data.values()(sample, source_probe) - mean
                : 0.0;
        }
    }
    Eigen::VectorXd sumsquares(probe_count);
    for (int probe = 0; probe < probe_count; probe++) {
        sumsquares(probe) = centered.row(probe).dot(centered.row(probe));
    }
    whole_correlation = centered * centered.transpose();
    for (int row = 0; row < probe_count; row++) {
        whole_correlation(row, row) = 1.0;
        for (int column = row + 1; column < probe_count; column++) {
            whole_correlation(row, column) = whole_correlation(column, row) =
                whole_correlation(row, column) / sqrt(sumsquares(row) * sumsquares(column));
        }
    }
}

void ECOJO::joint(const vector<int> &selected, Eigen::VectorXd &b_joint, Eigen::VectorXd &se_joint, Eigen::VectorXd &p_joint) {
    int count = selected.size();
    Eigen::VectorXd b(count), n(count);
    for (int index = 0; index < count; index++) {
        b(index) = b_values(selected[index]);
        n(index) = n_values(selected[index]);
    }
    b_joint = correlation * b;
    se_joint = Eigen::VectorXd::Constant(count, -1.0);
    p_joint = Eigen::VectorXd::Constant(count, 2.0);
    for (int index = 0; index < count; index++) {
        se_joint(index) = correlation(index, index) / n(index);
        if (se_joint(index) > 1e-30) {
            se_joint(index) = sqrt(se_joint(index));
            double z_value = b_joint(index) / se_joint(index);
            p_joint(index) = StatLib::pchisq(z_value * z_value, 1);
        }
    }
}

void ECOJO::conditional(const vector<int> &selected, const vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond) {
    Eigen::VectorXd b1(selected.size()), b2(remain.size());
    for (size_t index = 0; index < selected.size(); index++) b1(index) = b_values(selected[index]);
    for (size_t index = 0; index < remain.size(); index++) b2(index) = b_values(remain[index]);
    b_cond = Eigen::VectorXd::Constant(remain.size(), 0.0);
    se_cond = Eigen::VectorXd::Constant(remain.size(), -1.0);
    p_cond = Eigen::VectorXd::Constant(remain.size(), 2.0);
    for (size_t index = 0; index < remain.size(); index++) {
        Eigen::VectorXd adjusted = correlation * cross_correlation.col(index);
        b_cond(index) = b2(index) - adjusted.dot(b1);
        se_cond(index) = (1.0 - cross_correlation.col(index).dot(adjusted)) / n_values(remain[index]);
        if (se_cond(index) > 1e-30) {
            se_cond(index) = sqrt(se_cond(index));
            double z_value = b_cond(index) / se_cond(index);
            p_cond(index) = StatLib::pchisq(z_value * z_value, 1);
        }
    }
}

bool ECOJO::init_correlation(const vector<int> &selected) {
    correlation.resize(selected.size(), selected.size());
    for (size_t row = 0; row < selected.size(); row++) {
        for (size_t column = 0; column < selected.size(); column++) {
            correlation(row, column) = whole_correlation(selected[row], selected[column]);
        }
    }
    invert_correlation();
    return (1.0 - Eigen::VectorXd::Constant(selected.size(), 1.0).array() / correlation.diagonal().array()).maxCoeff() <= collinear;
}

void ECOJO::init_cross_correlation(const vector<int> &selected, const vector<int> &remain) {
    cross_correlation.resize(selected.size(), remain.size());
    for (size_t row = 0; row < selected.size(); row++) {
        for (size_t column = 0; column < remain.size(); column++) {
            cross_correlation(row, column) = whole_correlation(selected[row], remain[column]);
        }
    }
}

bool ECOJO::insert_correlation(const vector<int> &selected, int insert_index) {
    Eigen::MatrixXd previous = correlation;
    vector<int> indices = selected;
    indices.push_back(insert_index);
    std::stable_sort(indices.begin(), indices.end());
    correlation.resize(indices.size(), indices.size());
    for (size_t row = 0; row < indices.size(); row++) {
        for (size_t column = 0; column < indices.size(); column++) {
            correlation(row, column) = whole_correlation(indices[row], indices[column]);
        }
    }
    invert_correlation();
    if ((1.0 - Eigen::VectorXd::Constant(indices.size(), 1.0).array() / correlation.diagonal().array()).maxCoeff() > collinear) {
        correlation = previous;
        return false;
    }
    return true;
}

void ECOJO::erase_correlation(const vector<int> &selected) {
    correlation.resize(selected.size(), selected.size());
    for (size_t row = 0; row < selected.size(); row++) {
        for (size_t column = 0; column < selected.size(); column++) {
            correlation(row, column) = whole_correlation(selected[row], selected[column]);
        }
    }
    invert_correlation();
}

void ECOJO::invert_correlation() {
    double logdet = 0.0;
    if (!MatrixOperations::inverse_llt(correlation, logdet)) {
        LOGGER << "Note: no solution to LDLT decomposition. Switching to LU decomposition." << std::endl;
        correlation = correlation.lu().solve(Eigen::MatrixXd::Identity(correlation.rows(), correlation.cols()));
    }
}

void ECOJO::stepwise(vector<int> &selected, vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond) {
    vector<double> p_buffer(p_values.data(), p_values.data() + p_values.size());
    int best = std::min_element(p_buffer.begin(), p_buffer.end()) - p_buffer.begin();
    if (p_buffer[best] >= p_cutoff) return;
    selected.push_back(best);
    for (int index = 0; index < include.size(); index++) {
        if (index != best) remain.push_back(index);
    }
    int previous_count = 0;
    init_correlation(selected);
    init_cross_correlation(selected, remain);
    if (p_cutoff > 1e-3) LOGGER << "Performing forward model selection because the significance level is too low..." << std::endl;
    while (!remain.empty()) {
        if (select_entry(selected, remain, b_cond, se_cond, p_cond)) {
            if (p_cutoff <= 1e-3) select_stay(selected, b_cond, se_cond, p_cond);
            init_cross_correlation(selected, remain);
        } else {
            break;
        }
        if (selected.size() % 5 == 0 && selected.size() > previous_count) {
            LOGGER << selected.size() << " associated probes have been selected." << std::endl;
        }
        if (selected.size() > previous_count) previous_count = selected.size();
    }
    if (p_cutoff > 1e-3) {
        LOGGER << "Performing backward elimination..." << std::endl;
        select_stay(selected, b_cond, se_cond, p_cond);
    }
    LOGGER << "Finally, " << selected.size() << " associated probes are selected." << std::endl;
}

bool ECOJO::select_entry(vector<int> &selected, vector<int> &remain, Eigen::VectorXd &b_cond, Eigen::VectorXd &se_cond, Eigen::VectorXd &p_cond) {
    conditional(selected, remain, b_cond, se_cond, p_cond);
    vector<double> p_buffer(p_cond.data(), p_cond.data() + p_cond.size());
    while (!p_buffer.empty()) {
        int best = std::min_element(p_buffer.begin(), p_buffer.end()) - p_buffer.begin();
        if (p_buffer[best] >= p_cutoff) {
            init_cross_correlation(selected, remain);
            conditional(selected, remain, b_cond, se_cond, p_cond);
            return false;
        }
        if (insert_correlation(selected, remain[best])) {
            selected.push_back(remain[best]);
            std::stable_sort(selected.begin(), selected.end());
            remain.erase(remain.begin() + best);
            return true;
        }
        p_buffer.erase(p_buffer.begin() + best);
        remain.erase(remain.begin() + best);
    }
    return false;
}

void ECOJO::select_stay(vector<int> &selected, Eigen::VectorXd &b_joint, Eigen::VectorXd &se_joint, Eigen::VectorXd &p_joint) {
    while (!selected.empty()) {
        joint(selected, b_joint, se_joint, p_joint);
        vector<double> p_buffer(p_joint.data(), p_joint.data() + p_joint.size());
        int worst = std::max_element(p_buffer.begin(), p_buffer.end()) - p_buffer.begin();
        if (p_joint(worst) > p_cutoff) {
            selected.erase(selected.begin() + worst);
            erase_correlation(selected);
        } else {
            break;
        }
    }
}

void ECOJO::select(double p_cutoff_value, double collinear_value) {
    p_cutoff = p_cutoff_value;
    collinear = collinear_value;
    vector<int> selected, remain;
    Eigen::VectorXd b_cond, se_cond, p_cond;
    LOGGER << std::endl;
    LOGGER << "Performing stepwise model selection on " << include.size() << " probes to select association signals ... (p-value cutoff = " << p_cutoff << "; ";
    LOGGER << "collinearity cutoff = " << collinear << ")" << std::endl;
    stepwise(selected, remain, b_cond, se_cond, p_cond);
    if (selected.empty()) {
        LOGGER << "No probe has been selected." << std::endl;
        return;
    }
    Eigen::VectorXd b_joint, se_joint, p_joint;
    LOGGER << "Performing joint analysis on all the " << selected.size() << " selected signals ..." << std::endl;
    if (selected.size() >= sample_count) {
        LOGGER.e(0, "too many probes. The number of probes in a joint analysis should not be larger than the sample size.");
    }
    joint(selected, b_joint, se_joint, p_joint);
    write_selection(selected, b_joint, se_joint, p_joint);
}

void ECOJO::write_selection(const vector<int> &selected, const Eigen::VectorXd &b_joint, const Eigen::VectorXd &se_joint, const Eigen::VectorXd &p_joint) {
    string filename = output + ".slct.ecojo";
    LOGGER << "Saving the " << selected.size() << " independent signals to [" + filename + "] ..." << std::endl;
    std::ofstream output_file(filename);
    if (!output_file) LOGGER.e(0, "cannot open file [" + filename + "] to write.");
    output_file << "Probe\tb\tse\tz\tn\tbJ\tbJ_se\tzJ\tpJ" << std::endl;
    for (size_t index = 0; index < selected.size(); index++) {
        int probe = selected[index];
        output_file << probe_names[include[probe]] << "\t" << b_values(probe) << "\t" << se_values(probe)
                    << "\t" << z_values(probe) << "\t" << n_values(probe) << "\t" << b_joint(index)
                    << "\t" << se_joint(index) << "\t" << b_joint(index) / se_joint(index) << "\t" << p_joint(index) << "\t" << std::endl;
    }
}

void ECOJO::blup(double lambda) {
    LOGGER << "\nPerforming joint analysis on all the " << include.size() << " probes ..." << std::endl;
    double mean_n = n_values.mean();
    double diagonal = 1.0 + include.size() * (1.0 / lambda - 1.0) / mean_n;
    for (int index = 0; index < include.size(); index++) whole_correlation(index, index) = diagonal;
    double logdet = 0.0;
    if (!MatrixOperations::inverse_llt(whole_correlation, logdet)) {
        LOGGER << "Note: no solution to LDLT decomposition. Switching to LU decomposition." << std::endl;
        whole_correlation = whole_correlation.lu().solve(Eigen::MatrixXd::Identity(include.size(), include.size()));
    }
    Eigen::VectorXd b_blup = whole_correlation * b_values;
    string filename = output + ".blup.ecojo";
    LOGGER << "Saving the BLUP analysis result of " << include.size() << " probes to [" + filename + "] ..." << std::endl;
    std::ofstream output_file(filename);
    if (!output_file) LOGGER.e(0, "cannot open file [" + filename + "] to write.");
    output_file << "Probe\tb\tse\tz\tn\tb_blup" << std::endl;
    for (int index = 0; index < include.size(); index++) {
        output_file << probe_names[include[index]] << "\t" << b_values(index) << "\t" << se_values(index)
                    << "\t" << z_values(index) << "\t" << n_values(index) << "\t" << b_blup(index) << std::endl;
    }
}

int ECOJO::registerOption(map<string, vector<string>> &options_in) {
    options["out"] = options_in["out"][0];
    options["p"] = "5e-6";
    options["collinear"] = "0.9";
    action.clear();
    bool has_main = false;
    if (options_in.find("--ecojo-slct") != options_in.end()) {
        action = "select";
        has_main = true;
    }
    if (options_in.find("--ecojo-blup") != options_in.end()) {
        if (action.empty()) action = "blup";
        if (options_in["--ecojo-blup"].size() == 1) options["lambda"] = options_in["--ecojo-blup"][0];
        has_main = true;
    }
    if (!has_main) return 0;

    if (options_in.find("--efile") != options_in.end()) options["efile"] = options_in["--efile"][0];
    if (options_in.find("--e-cor") != options_in.end()) options["e-cor"] = options_in["--e-cor"][0];
    if (options_in.find("--ecojo") != options_in.end()) options["summary"] = options_in["--ecojo"][0];
    if (options_in.find("--ecojo-p") != options_in.end() && options_in["--ecojo-p"].size() == 1) options["p"] = options_in["--ecojo-p"][0];
    if (options_in.find("--ecojo-collinear") != options_in.end() && options_in["--ecojo-collinear"].size() == 1) {
        options["collinear"] = options_in["--ecojo-collinear"][0];
    }
    for (const string &flag : {string("--efile"), string("--e-cor"), string("--ecojo"), string("--ecojo-slct"), string("--ecojo-blup"), string("--ecojo-p"), string("--ecojo-collinear")}) {
        options_in.erase(flag);
    }
    return 1;
}

void ECOJO::processMain() {
    double p_cutoff = atof(options["p"].c_str());
    if (p_cutoff > 0.05 || p_cutoff <= 0.0) LOGGER.e(0, "\n  --ecojo-p should be within the range from 0 to 0.05.\n");
    double collinear = atof(options["collinear"].c_str());
    if (collinear > 1.0 || collinear < 0.01) LOGGER.e(0, "\n --ecojo-collinear should be within the range from 0.01 to 0.99.\n");
    if (options.find("summary") == options.end()) LOGGER.e(0, "please specify the EWAS summary data with --ecojo.");

    ECOJO ecojo(options["out"]);
    if (options.find("efile") != options.end()) {
        ExpressionData data(options["efile"]);
        ecojo.load_expression(data);
        ecojo.read_summary(options["summary"]);
        if (action == "select") {
            ecojo.select(p_cutoff, collinear);
        } else {
            double lambda = atof(options["lambda"].c_str());
            if (lambda < 0.01 || lambda > 0.99) LOGGER.e(0, "\n --ecojo-blup should be within the range from 0.01 to 0.99.\n");
            ecojo.blup(lambda);
        }
        return;
    }
    if (options.find("e-cor") == options.end()) LOGGER.e(0, "please specify expression data with --efile or --e-cor.");
    ecojo.load_correlation(options["e-cor"]);
    ecojo.read_summary(options["summary"]);
    double lambda = atof(options["lambda"].c_str());
    if (lambda < 0.01 || lambda > 0.99) LOGGER.e(0, "\n --ecojo-blup should be within the range from 0.01 to 0.99.\n");
    ecojo.blup(lambda);
}
