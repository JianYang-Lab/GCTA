#include "data/expression/ExpressionData.h"

#include "infrastructure/Logger.h"
#include <fstream>
#include <map>
#include <sstream>

using std::string;
using std::vector;

ExpressionData::ExpressionData(const string &filename) {
    std::ifstream input(filename);
    if (!input.is_open()) LOGGER.e(0, "cannot open the file [" + filename + "] to read.");
    LOGGER << "Reading gene expression / methylation data from [" + filename + "] ..." << std::endl;

    string line;
    std::getline(input, line);
    std::stringstream header(line);
    string name;
    if (!(header >> name)) LOGGER.e(0, "there needs to be at least 3 columns in the file [" + filename + "].");
    if (!(header >> name)) LOGGER.e(0, "there needs to be at least 3 columns in the file [" + filename + "].");
    std::map<string, int> probe_map;
    while (header >> name) {
        probe_map.insert(std::make_pair(name, static_cast<int>(probe_names.size())));
        if (probe_map.size() != probe_names.size() + 1) {
            LOGGER.e(0, "duplicated probe names found: \"" + name + "\".");
        }
        probe_names.push_back(name);
    }
    if (probe_names.empty()) LOGGER.e(0, "there needs to be at least 3 columns in the file [" + filename + "].");

    vector<vector<double>> rows;
    vector<vector<bool>> row_valid;
    std::map<string, int> sample_map;
    int line_number = 1;
    while (std::getline(input, line)) {
        line_number++;
        std::stringstream stream(line);
        string family_id, individual_id, value;
        if (!(stream >> family_id)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
        if (!(stream >> individual_id)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
        fid.push_back(family_id);
        pid.push_back(individual_id);
        string sample_id = family_id + ":" + individual_id;
        sample_map.insert(std::make_pair(sample_id, static_cast<int>(fid.size() - 1)));
        if (sample_map.size() != fid.size()) {
            LOGGER.e(0, "Duplicate individual ID found: \"" + family_id + "\t" + individual_id + "\".");
        }

        vector<double> row(probe_names.size());
        vector<bool> valid(probe_names.size());
        for (size_t column = 0; column < probe_names.size(); column++) {
            if (!(stream >> value)) LOGGER.e(0, "in line " + std::to_string(line_number) + ".");
            valid[column] = value != "-9";
            row[column] = valid[column] ? atof(value.c_str()) : 1e10;
        }
        rows.push_back(row);
        row_valid.push_back(valid);
    }

    probe_data.resize(rows.size(), probe_names.size());
    valid_data = row_valid;
    for (size_t row = 0; row < rows.size(); row++) {
        for (size_t column = 0; column < probe_names.size(); column++) {
            probe_data(row, column) = rows[row][column];
        }
    }
    LOGGER << "Expression data for " << probe_names.size() << " probes of " << rows.size()
           << " individuals have been included from the file [" + filename + "]." << std::endl;
}
