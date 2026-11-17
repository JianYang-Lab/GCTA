#ifndef GCTA2_NUMERICBACKEND_H
#define GCTA2_NUMERICBACKEND_H

namespace NumericBackend {

using DsyrkFunc = void (*)(const char *, const char *, const int *, const int *,
    const double *, const double *, const int *, const double *, double *, const int *);
using DgemmFunc = void (*)(const char *, const char *, const int *, const int *, const int *,
    const double *, const double *, const int *, const double *, const int *,
    const double *, double *, const int *);

extern DsyrkFunc dsyrk;
extern DgemmFunc dgemm;

void initialize();

}

#endif
