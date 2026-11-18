/*
 *  GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for MKL
 *
 * 2014 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"

/////////////////
// data functions

void gcta::make_XMat_mkl(float *X, bool grm_d_flag)
{
    if (_mu.empty()) calcu_mu();

    LOGGER << "Recoding genotypes (individual major mode) ..." << endl;
    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();

    if (!grm_d_flag) {
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            if (_dosage_flag) {
                for (j = 0; j < m; j++) {
                    if (_geno_dose[_keep[i]][_include[j]] < 1e5) {
                        if (_allele1[_include[j]] == _ref_A[_include[j]]) X[i * m + j] = _geno_dose[_keep[i]][_include[j]];
                        else X[i * m + j] = 2.0 - _geno_dose[_keep[i]][_include[j]];
                    } else X[i * m + j] = 1e6;
                }
                _geno_dose[i].clear();
            } else {
                for (j = 0; j < _include.size(); j++) {
                    if (!_snp_1[_include[j]][_keep[i]] || _snp_2[_include[j]][_keep[i]]) {
                        if (_allele1[_include[j]] == _ref_A[_include[j]]) X[i * m + j] = _snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]];
                        else X[i * m + j] = 2.0 - (_snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]]);
                    } else X[i * m + j] = 1e6;
                }
            }
        }
    } 
    else {
        #pragma omp parallel for private(j, k)
        for (i = 0; i < n; i++) {
            if (_dosage_flag) {
                for (j = 0; j < m; j++) {
                    if (_geno_dose[_keep[i]][_include[j]] < 1e5) {
                        k = i * m + j;
                        if (_allele1[_include[j]] == _ref_A[_include[j]]) X[k] = _geno_dose[_keep[i]][_include[j]];
                        else X[k] = 2.0 - _geno_dose[_keep[i]][_include[j]];
                        if (X[k] < 0.5) X[k] = 0.0;
                        else if (X[k] < 1.5) X[k] = _mu[_include[j]];
                        else X[k] = (2.0 * _mu[_include[j]] - 2.0);
                    } else X[k] = 1e6;
                }
                _geno_dose[i].clear();
            } 
            else {
                for (j = 0; j < _include.size(); j++) {
                    if (!_snp_1[_include[j]][_keep[i]] || _snp_2[_include[j]][_keep[i]]) {
                        k = i * m + j;
                        if (_allele1[_include[j]] == _ref_A[_include[j]]) X[k] = _snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]];
                        else X[k] = 2.0 - (_snp_1[_include[j]][_keep[i]] + _snp_2[_include[j]][_keep[i]]);
                        if (X[k] < 0.5) X[k] = 0.0;
                        else if (X[k] < 1.5) X[k] = _mu[_include[j]];
                        else X[k] = (2.0 * _mu[_include[j]] - 2.0);
                    } else X[i * m + j] = 1e6;
                }
            }
        }
    }
}

void gcta::std_XMat_mkl(float *X, vector<double> &sd_SNP, bool grm_xchr_flag, bool miss_with_mu, bool divid_by_std) {
    if (_mu.empty()) calcu_mu();

    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();
    sd_SNP.clear();
    sd_SNP.resize(m);
    if (_dosage_flag) {
#pragma omp parallel for private(i)
        for (j = 0; j < m; j++) {
            for (i = 0; i < n; i++) {
                double d_buf = X[i * m + j] - _mu[_include[j]];
                sd_SNP[j] += d_buf*d_buf;
            }
            sd_SNP[j] /= (n - 1.0);
        }
    } else {
        for (j = 0; j < m; j++) sd_SNP[j] = _mu[_include[j]]*(1.0 - 0.5 * _mu[_include[j]]);
    }
    if (divid_by_std) {
        for (j = 0; j < m; j++) {
            if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
            else sd_SNP[j] = sqrt(1.0 / sd_SNP[j]);
        }
    }

#pragma omp parallel for private(j, k)
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            k = i * m + j;
            if (X[k] < 1e5) {
                X[k] -= _mu[_include[j]];
                if (divid_by_std) X[k] *= sd_SNP[j];
            } else if (miss_with_mu) X[k] = 0.0;
        }
    }

    if (!grm_xchr_flag) return;
    // for the X-chromosome
    check_sex();
    double f_buf = sqrt(0.5);

#pragma omp parallel for private(j, k)
    for (i = 0; i < n; i++) {
        if (_sex[_keep[i]] == 1) {
            for (j = 0; j < m; j++) {
                k = i * m + j;
                if (X[k] < 1e5) X[k] *= f_buf;
                else if (miss_with_mu) X[k] = 0.0;
            }
        }
    }
}

void gcta::std_XMat_d_mkl(float *X, vector<double> &sd_SNP, bool miss_with_mu, bool divid_by_std) {
    if (_mu.empty()) calcu_mu();

    unsigned long i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();
    sd_SNP.clear();
    sd_SNP.resize(m);
    if (_dosage_flag) {
#pragma omp parallel for private(i)
        for (j = 0; j < m; j++) {
            for (i = 0; i < n; i++) {
                double d_buf = (X[i * m + j] - _mu[_include[j]]);
                sd_SNP[j] += d_buf*d_buf;
            }
            sd_SNP[j] /= (n - 1.0);
        }
    } else {
        for (j = 0; j < m; j++) sd_SNP[j] = _mu[_include[j]]*(1.0 - 0.5 * _mu[_include[j]]);
    }
    if (divid_by_std) {
        for (j = 0; j < m; j++) {
            if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
            else sd_SNP[j] = 1.0 / sd_SNP[j];
        }
    } else {
        for (j = 0; j < m; j++) sd_SNP[j] = sd_SNP[j] * sd_SNP[j];
    }
    vector<double> psq(m);
    for (j = 0; j < m; j++) psq[j] = 0.5 * _mu[_include[j]] * _mu[_include[j]];

#pragma omp parallel for private(j, k)
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            k = i * m + j;
            if (X[k] < 1e5) {
                X[k] -= psq[j];
                if (divid_by_std) X[k] *= sd_SNP[j];
            } else if (miss_with_mu) X[k] = 0.0;
        }
    }
}

////////////
// GRM
void gcta::make_grm_mkl(bool grm_d_flag, bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, bool mlmassoc, bool diag_f3_flag)
{
    if (!grm_d_flag && grm_xchr_flag) check_chrX();
    else check_autosome();

    unsigned long i = 0, j = 0, k = 0, l = 0, n = _keep.size(), m = _include.size();
    _geno_mkl = new float[n * m]; // alloc memory to X matrix

    make_XMat_mkl(_geno_mkl, grm_d_flag);
    vector<double> sd_SNP, sd_SNP_buf;
    if (grm_mtd == 0) {
        if (grm_d_flag) std_XMat_d_mkl(_geno_mkl, sd_SNP, false, true);
        else std_XMat_mkl(_geno_mkl, sd_SNP, grm_xchr_flag, false, true);
    } 
    else {
        if (grm_d_flag) std_XMat_d_mkl(_geno_mkl, sd_SNP, false, false);
        else std_XMat_mkl(_geno_mkl, sd_SNP, grm_xchr_flag, false, false);
    }

    if (!mlmassoc) LOGGER << "\nCalculating the" << ((grm_d_flag) ? " dominance" : "") << " genetic relationship matrix (GRM)" << (grm_xchr_flag ? " for the X chromosome" : "") << (_dosage_flag ? " using imputed dosage data" : "") << " ... (Note: default speed-optimized mode, may use huge RAM)" << endl;
    else LOGGER << "\nCalculating the genetic relationship matrix (GRM) ... " << endl;

    // count the number of missing genotypes
    vector< vector<int> > miss_pos(n);
    bool * X_bool = new bool[n * m];
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            k = i * m + j;
            if (_geno_mkl[k] < 1e5) X_bool[k] = true;
            else {
                _geno_mkl[k] = 0.0;
                miss_pos[i].push_back(j);
                X_bool[k] = false;
            }
        }
    }

    // Calculate A_N matrix
    _grm_N.resize(n, n);
    #pragma omp parallel for private(j, k)
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            int miss_j = 0;
            for (k = 0; k < miss_pos[j].size(); k++) miss_j += (int) X_bool[i * m + miss_pos[j][k]];
            _grm_N(i,j) = m - miss_pos[i].size() - miss_j;
        }
    }

    // Calculate sum of LD weights
    if (grm_mtd == 1) {
        double denominator = 0.0;
        for (j = 0; j < m; j++) denominator += sd_SNP[j];
        denominator = denominator / (double)m; 
        #pragma omp parallel for private(j, k)
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) _grm_N(i,j) = _grm_N(i,j) * denominator;
        }    
    }
    
    // Calcuate WW'
    _grm_mkl = new float[n * n]; // alloc memory to A
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, m, 1.0, _geno_mkl, m, _geno_mkl, m, 0.0, _grm_mkl, n);

    // re-calcuate the diagonals (Fhat3+1)
    if (diag_f3_flag) {
        #pragma omp parallel for private(j,k,l)
        for (i = 0; i < n; i++) {
            l = i * n + i;
            _grm_mkl[l] = 0.0;
            for (j = 0; j < m; j++) {
                k = i * m + j;
                _grm_mkl[l] += _geno_mkl[k]*(_geno_mkl[k]+(_mu[_include[j]] - 1.0) * sd_SNP[j]);
            }
        }
    }

    // Calculate A matrix
    #pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        for (j = 0; j <= i; j++) {
            if (_grm_N(i,j) > 0.0) _grm_mkl[i * n + j] /= _grm_N(i,j);
            else _grm_mkl[i * n + j] = 0.0;
        }
    }

    if (inbred) {
        #pragma omp parallel for private(j)
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) _grm_mkl[i * n + j] *= 0.5;
        }
    }

    if (mlmassoc && grm_mtd == 0) {
        for (j = 0; j < m; j++) {
            if (fabs(sd_SNP[j]) < 1.0e-50) sd_SNP[j] = 0.0;
            else sd_SNP[j] = 1.0 / sd_SNP[j];
        }
        #pragma omp parallel for private(j, k)
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                k = i * m + j;
                if (_geno_mkl[k] < 1e5) _geno_mkl[k] *= sd_SNP[j];
                else _geno_mkl[k] = 0.0;
            }
        }
        delete[] X_bool;
    } else {
        // Output A_N and A
        string out_buf = _out;
        if (grm_d_flag) _out += ".d";
        output_grm_mkl(_grm_mkl, output_bin);
        _out = out_buf;

        // free memory
        delete[] _geno_mkl;
        delete[] X_bool;
        delete[] _grm_mkl;
    }
}

void gcta::output_grm_mkl(float* A, bool output_grm_bin)
{
    unsigned long i = 0, j = 0, n = _keep.size();
    string grm_file;

    if (output_grm_bin) {
        // Save matrix A in binary file
        grm_file = _out + ".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out | ios::binary);
        if (!A_Bin) LOGGER.e(0, "cannot open the file [" + grm_file + "] to write.");
        int size = sizeof (float);
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) A_Bin.write((char*) &(A[i * n + j]), size);
        }
        A_Bin.close();
        LOGGER << "GRM of " << n << " individuals has been saved in the file [" + grm_file + "] (in binary format)." << endl;

        string grm_N_file = _out + ".grm.N.bin";
        fstream N_Bin(grm_N_file.c_str(), ios::out | ios::binary);
        if (!N_Bin) LOGGER.e(0, "cannot open the file [" + grm_N_file + "] to write.");
        size = sizeof (float);
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) N_Bin.write((char*) &(_grm_N(i,j)), size);
        }
        N_Bin.close();
        LOGGER << "Number of SNPs to calculate the genetic relationship between each pair of individuals has been saved in the file [" + grm_N_file + "] (in binary format)." << endl;
    } else {
        // Save A matrix in txt format
        grm_file = _out + ".grm.gz";
        gzofstream zoutf;
        zoutf.open(grm_file.c_str());
        if (!zoutf.is_open()) LOGGER.e(0, "cannot open the file [" + grm_file + "] to write.");
        LOGGER << "Saving the genetic relationship matrix to the file [" + grm_file + "] (in compressed text format)." << endl;
        zoutf.setf(ios::scientific);
        zoutf.precision(6);
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) zoutf << i + 1 << '\t' << j + 1 << '\t' << _grm_N(i,j) << '\t' << A[i * n + j] << endl;
        }
        zoutf.close();
        LOGGER << "The genetic relationship matrix has been saved in the file [" + grm_file + "] (in compressed text format)." << endl;
    }

    string famfile = _out + ".grm.id";
    ofstream Fam(famfile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open the file [" + famfile + "] to write.");
    for (i = 0; i < n; i++) Fam << _fid[_keep[i]] + "\t" + _pid[_keep[i]] << endl;
    Fam.close();
    LOGGER << "IDs for the GRM file [" + grm_file + "] have been saved in the file [" + famfile + "]." << endl;
}

///////////
// reml

bool gcta::comput_inverse_logdet_LDLT_mkl(eigenMatrix &Vi, double &logdet)
{
    //LOGGER << "LDLT" << endl;
    uint64_t n = Vi.cols();
    double* Vi_mkl = new double[n * n];
    //float* Vi_mkl=new float[n*n];

    #pragma omp parallel for
    for (uint64_t j = 0; j < n; j++) {
        for (uint64_t i = 0; i < n; i++) {
            Vi_mkl[i * n + j] = Vi(i, j);
        }
    }
    //LOGGER << "Finished copy" << endl;

    // MKL's Cholesky decomposition
    int info = 0, int_n = (int) n;
    char uplo = 'L';
#if GCTA_CPU_x86
    dpotrf(&uplo, &int_n, Vi_mkl, &int_n, &info);
#else
    dpotrf_(&uplo, &int_n, Vi_mkl, &int_n, &info);
#endif
    //LOGGER << "Finished decompose" << endl;
    //spotrf( &uplo, &n, Vi_mkl, &n, &info );
    if (info < 0){
        LOGGER.e(0, "Cholesky decomposition failed. Invalid values found in the matrix.\n");
    }else if (info > 0){
        delete[] Vi_mkl;
        return false;
    }else {
        logdet = 0.0;
        for (uint64_t i = 0; i < n; i++) {
            double d_buf = Vi_mkl[i * n + i];
            logdet += log(d_buf * d_buf);
        }

        //LOGGER << "start inverse" << endl;
        // Calcualte V inverse
#if GCTA_CPU_x86
        dpotri(&uplo, &int_n, Vi_mkl, &int_n, &info);
#else
        dpotri_(&uplo, &int_n, Vi_mkl, &int_n, &info);
#endif
        //LOGGER << "Inverse finished" << endl;
        //spotri( &uplo, &n, Vi_mkl, &n, &info );
        if (info < 0){
            LOGGER.e(0, "invalid values found in the variance-covariance (V) matrix.\n");
        }else if (info > 0){
            delete[] Vi_mkl;
            return false;
        }else {
            #pragma omp parallel for
            for (uint64_t j = 0; j < n; j++) {
                for (uint64_t i = 0; i <= j; i++) Vi(i, j) = Vi(j, i) = Vi_mkl[i * n + j];
            }
            //LOGGER << "IDLT finished" << endl;
        }
    }

    // free memory
    delete[] Vi_mkl;

    return true;

}

bool gcta::comput_inverse_logdet_LU_mkl(eigenMatrix &Vi, double &logdet)
{
    unsigned long n = Vi.cols();
    double* Vi_mkl = new double[n * n];

    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            Vi_mkl[i * n + j] = Vi(i, j);
        }
    }

    int N = (int) n;
    int *IPIV = new int[n + 1];
    int LWORK = N*N;
    double *WORK = new double[n * n];
    int INFO;
#if GCTA_CPU_x86
    dgetrf(&N, &N, Vi_mkl, &N, IPIV, &INFO);
#else
    dgetrf_(&N, &N, Vi_mkl, &N, IPIV, &INFO);
#endif
    if (INFO < 0) LOGGER.e(0, "LU decomposition failed. Invalid values found in the matrix.\n");
    else if (INFO > 0) {
        delete[] Vi_mkl;
        delete[] IPIV;
        delete[] WORK;
        return false;
    } else {
        logdet = 0.0;
        for (uint64_t i = 0; i < n; i++) {
            double d_buf = Vi_mkl[i * n + i];
            logdet += log(fabs(d_buf));
        }

        // Calcualte V inverse
#if GCTA_CPU_x86
        dgetri(&N, Vi_mkl, &N, IPIV, WORK, &LWORK, &INFO);
#else
        dgetri_(&N, Vi_mkl, &N, IPIV, WORK, &LWORK, &INFO);
#endif
        if (INFO < 0){
            LOGGER.e(0, "invalid values found in the variance-covariance (V) matrix.\n");
        }else if (INFO > 0){
            delete[] Vi_mkl;
            delete[] IPIV;
            delete[] WORK;
            return false;
        }else {
            #pragma omp parallel for
            for (uint64_t j = 0; j < n; j++) {
                for (uint64_t i = 0; i <= j; i++) Vi(i, j) = Vi(j, i) = Vi_mkl[i * n + j];
            }
        }
    }

    // free memory
    delete[] Vi_mkl;
    delete[] IPIV;
    delete[] WORK;

    return true;
}

bool gcta::comput_inverse_logdet_LU_mkl_array(int n, float *Vi, double &logdet) {
    double* Vi_mkl = new double[n * n];

    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            Vi_mkl[i * n + j] = Vi[i * n + j];
        }
    }

    int N = (int) n;
    int *IPIV = new int[n + 1];
    int LWORK = N*N;
    double *WORK = new double[n * n];
    int INFO;
#if GCTA_CPU_x86
    dgetrf(&N, &N, Vi_mkl, &N, IPIV, &INFO);
#else
    dgetrf_(&N, &N, Vi_mkl, &N, IPIV, &INFO);    
#endif
    if (INFO < 0) LOGGER.e(0, "LU decomposition failed. Invalid values found in the matrix.\n");
    else if (INFO > 0) {
        // free memory
        delete[] Vi_mkl;
        delete[] IPIV;
        delete[] WORK;

        return (false); //Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
    }else {
        logdet = 0.0;
        for (uint64_t i = 0; i < n; i++) {
            double d_buf = Vi_mkl[i * n + i];
            logdet += log(fabs(d_buf));
        }

        // Calcualte V inverse
#if GCTA_CPU_x86
        dgetri(&N, Vi_mkl, &N, IPIV, WORK, &LWORK, &INFO);
#else
        dgetri_(&N, Vi_mkl, &N, IPIV, WORK, &LWORK, &INFO);
#endif
        if (INFO < 0) LOGGER.e(0, "invalid values found in the variance-covariance (V) matrix.\n");
        else if (INFO > 0) {
            // free memory
            delete[] Vi_mkl;
            delete[] IPIV;
            delete[] WORK;
            return (false); // Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
        }else {
            #pragma omp parallel for
            for (uint64_t j = 0; j < n; j++) {
                for (uint64_t i = 0; i < n; i++) Vi[i * n + j] = Vi_mkl[i * n + j];
            }
        }
    }

    // free memory
    delete[] Vi_mkl;
    delete[] IPIV;
    delete[] WORK;

    return true;
}

//////////////
// LD

void gcta::LD_pruning_mkl(double rsq_cutoff, int wind_size) {
    check_autosome();

    unsigned long i = 0, j = 0, k = 0, l = 0, n = _keep.size(), m = _include.size();
    _geno_mkl = new float[n * m]; // alloc memory to X matrix
    make_XMat_mkl(_geno_mkl, false);
    vector<double> sd_SNP;
    std_XMat_mkl(_geno_mkl, sd_SNP, false, true, true);

    LOGGER << "\nPruning SNPs for LD ..." << endl;
    vector<int> brk_pnt1, brk_pnt2, brk_pnt3;
    get_ld_blk_pnt(brk_pnt1, brk_pnt2, brk_pnt3, wind_size*2);

    vector<int> rm_snp_indx;
    LD_pruning_blk_mkl(_geno_mkl, brk_pnt1, rsq_cutoff, rm_snp_indx);
    if (brk_pnt2.size() > 1) LD_pruning_blk_mkl(_geno_mkl, brk_pnt2, rsq_cutoff, rm_snp_indx);
    stable_sort(rm_snp_indx.begin(), rm_snp_indx.end());
    rm_snp_indx.erase(unique(rm_snp_indx.begin(), rm_snp_indx.end()), rm_snp_indx.end());
    int m_sub = rm_snp_indx.size();
    vector<string> rm_snp_name(m_sub);
#pragma omp parallel for
    for (i = 0; i < m_sub; i++) rm_snp_name[i] = _snp_name[_include[rm_snp_indx[i]]];
    update_id_map_rm(rm_snp_name, _snp_name_map, _include);
    m = _include.size();

    LOGGER << "After LD-pruning, " << m << " SNPs remain." << endl;
    string pruned_file = _out + ".prune.in";
    ofstream oprune(pruned_file.data());
    for (i = 0; i < m; i++) oprune << _snp_name[_include[i]] << endl;
    oprune << endl;
    LOGGER << "The list of " << m << " LD-pruned SNPs (pruned in) has been saved in the file [" + pruned_file + "]." << endl;
}

void gcta::LD_pruning_blk_mkl(float *X, vector<int> &brk_pnt, double rsq_cutoff, vector<int> &rm_snp_ID1)
{
    unsigned long i = 0, j = 0, k = 0, l = 0, n = _keep.size(), m = _include.size(), size = 0;

    for (i = 0; i < brk_pnt.size() - 1; i++) {
        if (_chr[_include[brk_pnt[i]]] != _chr[_include[brk_pnt[i + 1]]]) continue;
        size = brk_pnt[i + 1] - brk_pnt[i] + 1;
        if (size < 3) continue;

        float *X_sub = new float[n * size];
        #pragma omp parallel for private(k, l)
        for (j = 0; j < n; j++) {
            for (k = 0, l = brk_pnt[i]; k < size; k++, l++) X_sub[j * size + k] = X[j * m + l];
        }

        float *rsq_sub = new float[size * size];
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, size, size, n, 1.0, X_sub, size, X_sub, size, 0.0, rsq_sub, size);
        #pragma omp parallel for private(k, l)
        for (j = 0; j < size; j++) {
            for (k = 0; k < size; k++) {
                l = j * size + k;
                rsq_sub[l] /= (float) n;
                rsq_sub[l] = rsq_sub[l] * rsq_sub[l];
            }
        }
        vector<int> rm_snp_buf;
        rm_cor_snp((int) size, brk_pnt[i], rsq_sub, rsq_cutoff, rm_snp_buf);
        rm_snp_ID1.insert(rm_snp_ID1.end(), rm_snp_buf.begin(), rm_snp_buf.end());

        delete[] X_sub;
        delete[] rsq_sub;
    }
}
