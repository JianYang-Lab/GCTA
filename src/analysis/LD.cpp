#include "analysis/LD.h"

#include "analysis/StatLib.h"

#include "infra/constants.hpp"
#include "infra/cpu_f77blas.h"
#include "infra/Logger.h"
#include "infra/utils.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdio>
#include <functional>
#include <numeric>
#include <string>

map<string, string> LD::options;
map<string, int> LD::options_i;
vector<string> LD::processFunctions;
bool LD::chr_ends = false;
unique_ptr<double[]> LD::geno_buffer[2];
int LD::cur_buffer = 0;
uint64_t LD::cur_buffer_offset[2] = {0, 0};
uint32_t LD::num_indi = 0;

LD::LD(Geno * geno){
    this->geno = geno;
    num_indi = geno->pheno->count_keep();
    ld_window = options_i["LD_window"] * 1000;
    is_r2 = false;
    cur_process_marker_index = 0;
    if(options["method"] == "r2"){
        is_r2 = true;
    }

    h_ld = fopen(options["out"].c_str(), "wb");
    if(!h_ld){
        LOGGER.e(0, "can't open " + options["out"] + " for writing.");
    }
}

LD::~LD(){
    fclose(h_ld);
}

template <typename T>
struct static_cast_func
{
  template <typename T1> 
  T operator()(const T1& x) const { return static_cast<T>(x); }
};

void LD::calcLD(){
    // the current buffer
    int cacl_index_buffer;
    if(geno_buffer[!cur_buffer]){
        cacl_index_buffer = !cur_buffer;
    }else{
        cacl_index_buffer = cur_buffer;
    }
    char trans = 'T', notrans = 'N', uplo = 'L';
    double zero = 0.0;
    int nr = num_indi;
    int nc1 = cur_buffer_offset[cacl_index_buffer] / nr;
    double alpha = 1.0 / (nr - 1);
    double *ptr1 = geno_buffer[cacl_index_buffer].get();
    double *res1 = new double[nc1 * nc1];

    // GWAS LD is a band: only markers within the physical window are written.
    vector<uint32_t> band_sizes(nc1);
    for(int i = 0; i < nc1; i++){
        band_sizes[i] = geno->marker->getNextWindowSize(cur_process_marker_index + i, ld_window);
    }

#if GCTA_CPU_x86
    dsyrk(&uplo, &trans, &nc1, &nr, &alpha, ptr1, &nr, &zero, res1, &nc1);
#else
    dsyrk_(&uplo, &trans, &nc1, &nr, &alpha, ptr1, &nr, &zero, res1, &nc1);
#endif

    double *ptr2 = geno_buffer[!cacl_index_buffer].get();
    const int nc2 = ptr2 ? cur_buffer_offset[!cacl_index_buffer] / nr : 0;
    const int row_block = 64;
    vector<double> cross;
    vector<float> row;
    for(int r0 = 0; r0 < nc1; r0 += row_block){
        const int br = std::min(row_block, nc1 - r0);
        int need = 0;
        if(ptr2){
            for(int i = r0; i < r0 + br; i++){
                const uint32_t inside = nc1 - i;
                if(band_sizes[i] > inside){
                    need = std::max(need, static_cast<int>(std::min<uint32_t>(band_sizes[i] - inside, nc2)));
                }
            }
        }
        if(need){
            cross.resize(static_cast<size_t>(need) * br);
#if GCTA_CPU_x86
            dgemm(&trans, &notrans, &need, &br, &nr, &alpha, ptr2, &nr, ptr1 + (int64_t)r0 * nr, &nr, &zero, cross.data(), &need);
#else
            dgemm_(&trans, &notrans, &need, &br, &nr, &alpha, ptr2, &nr, ptr1 + (int64_t)r0 * nr, &nr, &zero, cross.data(), &need);
#endif
        }
        for(int j = 0; j < br; j++){
            const int i = r0 + j;
            const uint32_t inside = nc1 - i;
            const uint32_t cur_size = ptr2 ? band_sizes[i] : inside;
            if(cur_size < inside || cur_size - inside > static_cast<uint32_t>(need)){
                LOGGER.e(0, "inconsistent LD window boundaries.");
            }
            row.resize(cur_size);
            const double *within = res1 + static_cast<size_t>(i) * nc1 + i;
            std::transform(within, within + inside, row.data(), static_cast_func<float>());
            if(cur_size > inside){
                const double *outside = cross.data() + static_cast<size_t>(j) * need;
                std::transform(outside, outside + cur_size - inside, row.data() + inside, static_cast_func<float>());
            }
            if(fwrite(&cur_size, sizeof(uint32_t), 1, h_ld) != 1 ||
               fwrite(row.data(), sizeof(float), cur_size, h_ld) != cur_size){
                LOGGER.e(0, "can't write to " + options["out"] + ".");
            }
        }
    }
    cur_process_marker_index += nc1;
    delete[] res1;
    geno_buffer[cacl_index_buffer].reset(nullptr);
    fflush(h_ld);
}

void LD::readGeno(uint64_t *buf, int num_marker){
    double * ptr = geno_buffer[cur_buffer].get();
    uint64_t cur_offset = cur_buffer_offset[cur_buffer];
    ptr += cur_offset;

    #pragma omp parallel for schedule(dynamic) 
    for(int i = 0; i < num_marker; i++){
        geno->makeMarkerX(buf, i, ptr + i * num_indi, true, true);
    }
    cur_buffer_offset[cur_buffer] += (uint64_t) num_marker * num_indi;
}


int LD::registerOption(map<string, vector<string>>& options_in){
    int ret_val = 0;
    options["out"] = options_in["out"][0] + ".zld";

    string curFlag = "--ld-matrix";
    if(options_in.find(curFlag) != options_in.end()){
        processFunctions.push_back("ld_matrix");
        ret_val++;
        options_in.erase(curFlag);
    }

    // get the calculation methods, defalut LD r;
    options["method"] = "r";
    curFlag = "--r";
    if(options_in.find(curFlag) != options_in.end()){
        options["method"] = "r";
        options_in.erase(curFlag);
    }

    curFlag = "--r2";
    if(options_in.find(curFlag) != options_in.end()){
        options["method"] = "r2";
        options_in.erase(curFlag);
    }

    // set the LD window size Kb
    options_i["LD_window"] = 1000;
    curFlag = "--ld-wind";
    if(options_in.find(curFlag) != options_in.end()){
        if(options_in[curFlag].size() == 1){
            try{
                options_i["LD_window"] = std::stoi(options_in[curFlag][0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "LD window is not an integer.");
            }
        }
        options_in.erase(curFlag);
    }


    return ret_val;
}


void LD::processMain(){
    vector<function<void (uint64_t *, int)>> callBacks;
    for(auto &process_function : processFunctions){
        if(process_function == "ld_matrix"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.init_keep();
            geno.AFA1.resize(marker.count_extract());
            geno.countMarkers.resize(marker.count_extract());
            geno.init_AsyncBuffer();
            if(!geno.filterMAF()){
                vector<function<void (uint64_t *, int)>> freqCallBacks;
                freqCallBacks.push_back(std::bind(&Geno::freq64, &geno, _1, _2));
                geno.loop_64block(marker.get_extract_index(), freqCallBacks);
                geno.resetLoop();
                geno.init_AsyncBuffer();
            }
 
            LD ld(&geno);

            LOGGER.i(0, "Generating LD matrix...");
            uint32_t window = options_i["LD_window"] * 1000;
            uint32_t total_num_marker = marker.count_extract(); 
            uint32_t cur_index_marker = 0;
           callBacks.push_back(std::bind(&LD::readGeno, &ld, _1, _2));

            while(cur_index_marker < total_num_marker){
                cur_buffer = !cur_buffer;
                bool isX;
                vector<uint32_t> indices1 = marker.getNextWindowIndex(cur_index_marker, window, chr_ends, isX);
                geno_buffer[cur_buffer].reset(new double[indices1.size() * num_indi]);
                cur_buffer_offset[cur_buffer] = 0;
                geno.num_finished_markers = cur_index_marker;
                geno.init_AsyncBuffer();
                geno.loop_64block(indices1, callBacks, false);
                cur_index_marker += indices1.size();

                // if another buffer is NA,  not chr ends and still in range;
                if((!geno_buffer[!cur_buffer]) && (!chr_ends)){
                    continue;
                }
                //ld.deduceLD();
                ld.calcLD();
                LOGGER.p(0, to_string_precision(cur_index_marker * 100.0 / total_num_marker, 2) + "% finished.");
            }
            //last block if not chr ends;
            if(geno_buffer[cur_buffer]){
                ld.calcLD();
            }
        }
    }
}
