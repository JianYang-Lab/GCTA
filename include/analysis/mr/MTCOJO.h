#ifndef GCTA2_MTCOJO_H
#define GCTA2_MTCOJO_H

#include "analysis/mr/MRAnalysis.h"

class MTCOJO : public MRAnalysis {
public:
    explicit MTCOJO(std::string out) : MRAnalysis(std::move(out)) {}
    static int registerOption(std::map<std::string, std::vector<std::string>> &options_in);
    static void processMain();
};

#endif
