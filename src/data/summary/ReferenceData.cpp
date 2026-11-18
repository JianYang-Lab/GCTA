#include "data/summary/ReferenceData.h"

#include "infra/Logger.h"

#include <sstream>

using std::string;
using std::vector;

void ReferenceData::calcu_mu() {
    if (!loaded) load();
}

bool ReferenceData::make_XMat_subset(Eigen::MatrixXf &matrix, std::vector<int> &snp_indices, bool divide_by_std) {
    if (snp_indices.empty()) return false;
    calcu_mu();
    matrix.resize(_keep.size(), snp_indices.size());
    for (size_t column = 0; column < snp_indices.size(); column++) {
        int raw_index = _include[snp_indices[column]];
        const GenoBufItem &item = geno_items.at(raw_index);
        for (size_t sample = 0; sample < _keep.size(); sample++) {
            float value = item.geno[sample];
            if (!divide_by_std && item.sd > 1e-50) value /= sqrt(1.0 / item.sd);
            matrix(sample, column) = value;
        }
    }
    return true;
}

void ReferenceData::makex_eigenVector(int index, Eigen::VectorXd &values, bool resize, bool minus_2p) {
    calcu_mu();
    if (resize) values.resize(_keep.size());
    int raw_index = _include[index];
    const GenoBufItem &item = geno_items.at(raw_index);
    for (size_t sample = 0; sample < _keep.size(); sample++) {
        values(sample) = item.geno[sample];
        if (item.sd > 1e-50) values(sample) /= sqrt(1.0 / item.sd);
        if (!minus_2p) values(sample) += _mu[raw_index];
    }
}

V2ReferenceData::V2ReferenceData() {
    pheno = std::make_unique<Pheno>();
    marker = std::make_unique<Marker>();
    geno = std::make_unique<Geno>(pheno.get(), marker.get());
    geno->init_keep();

    int raw_count = marker->count_raw();
    _keep.resize(pheno->count_keep());
    for (int index = 0; index < pheno->count_keep(); index++) _keep[index] = index;
    for (int raw_index = 0; raw_index < raw_count; raw_index++) {
        string info = marker->get_marker(raw_index);
        std::stringstream stream(info);
        string chr, name, position, a1, a2;
        stream >> chr >> name >> position >> a1 >> a2;
        _chr.push_back(std::stoi(chr));
        _snp_name.push_back(name);
        _bp.push_back(std::stoi(position));
        _allele1.push_back(a1);
        _allele2.push_back(a2);
        _ref_A.push_back(a1);
        _snp_name_map.insert(std::make_pair(name, raw_index));
        _include.push_back(raw_index);
    }
}

void V2ReferenceData::load() {
    if (loaded) return;
    if (_include.empty()) LOGGER.e(0, "no SNP is retained for analysis.");
    vector<uint32_t> raw_include(_include.begin(), _include.end());
    marker->keep_raw_index(raw_include);
    vector<uint32_t> extracted(_include.size());
    for (size_t index = 0; index < extracted.size(); index++) extracted[index] = index;
    _mu.assign(marker->count_raw(), 0.0);
    geno_items.clear();
    vector<GenoBufItem> items;
    items.reserve(_include.size());
    auto callback = [this, &raw_include](uintptr_t *buffer, const vector<uint32_t> &indices) {
        for (size_t index = 0; index < indices.size(); index++) {
            GenoBufItem item;
            item.extractedMarkerIndex = indices[index];
            geno->getGenoDouble(buffer, index, &item);
            int raw_index = raw_include[indices[index]];
            _mu[raw_index] = 2.0 * item.af;
            geno_items[raw_index] = std::move(item);
        }
    };
    geno->loopDouble(extracted, 128, true, true, true, true, {callback}, false);
    loaded = true;
}
