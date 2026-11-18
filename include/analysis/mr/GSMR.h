#ifndef GCTA2_GSMR_H
#define GCTA2_GSMR_H

#include "analysis/mr/MRAnalysis.h"

class GSMR : public MRAnalysis {
public:
    explicit GSMR(std::string out) : MRAnalysis(std::move(out)) {}
    static int registerOption(std::map<std::string, std::vector<std::string>> &options_in);
    static void processMain();
};

#endif
