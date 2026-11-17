#include "NumericBackend.h"
#include "cpu_f77blas.h"
#include "Logger.h"
#include <cstring>

using std::string;

#if defined(__linux__) && GCTA_CPU_x86
#include <cpuid.h>
#include <dlfcn.h>
#endif

namespace NumericBackend {

#if GCTA_CPU_x86
DsyrkFunc dsyrk = &::dsyrk;
DgemmFunc dgemm = &::dgemm;
static std::string current_name = "MKL";
#else
DsyrkFunc dsyrk = &::dsyrk_;
DgemmFunc dgemm = &::dgemm_;
static std::string current_name = "SYSTEM";
#endif

#if defined(__linux__) && GCTA_CPU_x86
static void *openblas_handle = NULL;

static bool is_amd(){
    unsigned int eax, ebx, ecx, edx;
    if(!__get_cpuid(0, &eax, &ebx, &ecx, &edx)) return false;
    char vendor[13] = {};
    memcpy(vendor, &ebx, 4);
    memcpy(vendor + 4, &edx, 4);
    memcpy(vendor + 8, &ecx, 4);
    string vendor_string(vendor);
    return vendor_string == "AuthenticAMD" || vendor_string == "HygonGenuine";
}

static bool load_openblas(){
    if(openblas_handle) return true;
    void *handle = dlopen("libopenblaso.so.0", RTLD_NOW | RTLD_LOCAL);
    if(!handle) handle = dlopen("libopenblaso.so", RTLD_NOW | RTLD_LOCAL);
    if(!handle) return false;
    DsyrkFunc openblas_dsyrk = reinterpret_cast<DsyrkFunc>(dlsym(handle, "dsyrk_"));
    DgemmFunc openblas_dgemm = reinterpret_cast<DgemmFunc>(dlsym(handle, "dgemm_"));
    if(!openblas_dsyrk || !openblas_dgemm){
        dlclose(handle);
        return false;
    }
    openblas_handle = handle;
    dsyrk = openblas_dsyrk;
    dgemm = openblas_dgemm;
    current_name = "OPENBLAS";
    return true;
}
#endif

void initialize(){
    static bool initialized = false;
    if(initialized) return;
    initialized = true;
#if defined(__linux__) && GCTA_CPU_x86
    const char *requested = getenv("GCTA_NUMERIC_BACKEND");
    string choice = requested ? requested : "auto";
    if(choice != "auto" && choice != "mkl" && choice != "openblas"){
        LOGGER.e(0, "unknown numerical backend: " + choice + ". Expected auto, mkl, or openblas.");
    }
    bool use_openblas = choice == "openblas";
    if(choice == "auto" && is_amd()){
        LOGGER.i(0, "AMD CPU detected; keeping MKL as the conservative GRM default.");
    }
    if(use_openblas && !load_openblas()){
        LOGGER.e(0, "OpenBLAS was requested but libopenblaso with dsyrk/dgemm is unavailable.");
    }
#endif
    LOGGER.i(0, "GRM BLAS backend: " + current_name);
}

}
