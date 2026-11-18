#include "gcta.h"

#include "infra/Logger.h"

#include <algorithm>
#include <iterator>

const double eps = 1e-6;

void gcta::init_meta_snp_map(vector<string> snplist, map<string, int> &snp_name_map, vector<string> &snp_name, vector<int> &remain_snp) {
    int i=0, size=0, nsnp = snplist.size();

    snp_name_map.clear(); remain_snp.clear();
    remain_snp.resize(nsnp); 
    for(i=0; i<nsnp; i++) {
        snp_name_map.insert( pair<string,int>(snplist[i], i) );
        if(size==snp_name_map.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[i] + ".");
        size = snp_name_map.size();
        remain_snp[i] = i;
    }

    snp_name.clear(); snp_name = snplist;
}

void gcta::init_gwas_variable(vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq, eigenMatrix &snp_b, eigenMatrix &snp_se, eigenMatrix &snp_pval, eigenMatrix &n, int npheno, int nsnp) {
    int i = 0;
    // two alleles
    snp_a1.clear(); snp_a2.clear();
    snp_a1.resize(npheno); snp_a2.resize(npheno);
    for( i=0; i<npheno; i++) {
        snp_a1[i].resize(nsnp); snp_a2[i].resize(nsnp);
    }
    // allele frequency
    snp_freq.resize(nsnp, npheno);
    // effect size, se, p-value and sample size
    snp_b.resize(nsnp, npheno); snp_se.resize(nsnp, npheno);
    snp_pval.resize(nsnp, npheno); n.resize(nsnp, npheno);
}

void gcta::update_meta_snp_list(vector<string> &snplist, map<string, int> snp_id_map) {
    int i = 0, nsnpbuf = 0;
    vector<string> snpbuf(snplist);
    nsnpbuf = snplist.size();

    snplist.clear();
    for(i=0; i<nsnpbuf; i++) {
        if(snp_id_map.find(snpbuf[i]) == snp_id_map.end()) continue;
        snplist.push_back(snpbuf[i]);
    }
}

void gcta::update_meta_snp_map(vector<string> snplist, map<string, int> &snp_id_map, vector<string> &snp_id, vector<int> &snp_indx, bool indx_flag) {
    int i = 0, j = 0, size = 0, nsnp = snplist.size();
    map<string,int> snp_add_map;
    
    for(i=0, j=snp_id_map.size(); i<nsnp; i++) {
        // Check duplicated SNP ID
        snp_add_map.insert(pair<string,int>(snplist[i], i));
        if(size==snp_add_map.size()) LOGGER.e(0, "Duplicated SNP ID found: " + snplist[i] + ".");
        size = snp_add_map.size();

        if(snp_id_map.find(snplist[i]) != snp_id_map.end()) continue;
        snp_id_map.insert(pair<string, int>(snplist[i], j));
        snp_id.push_back(snplist[i]);
        if(indx_flag) snp_indx.push_back(j);
        j++;
    }
}

void gcta::update_meta_snp(map<string,int> &snp_name_map, vector<string> &snp_name, vector<int> &snp_remain) {
    int i = 0, nsnp = snp_remain.size();
    vector<int> snp_remain_buf;
    vector<string> snp_name_buf;

    snp_name_buf = snp_name;
    snp_remain_buf = snp_remain;
    // reset snp_name and snp_remain
    snp_name.clear(); snp_name.resize(nsnp);
    snp_remain.clear(); snp_remain.resize(nsnp);
    for(i = 0; i<nsnp; i++) {
        snp_name[i] = snp_name_buf[snp_remain_buf[i]];
        snp_remain[i] = i;
    }
    // reset snp_name_map
    snp_name_map.clear();
    for(i=0; i<nsnp; i++) 
        snp_name_map.insert(pair<string, int>(snp_name[i], i));
}

void update_mtcojo_snp_kp(const vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &remain_snp_indx) {
    
    int i=0, nsnpbuf = adjsnps.size();
    map<string,int> snp_id_map_buf(snp_id_map);
    for(i=0; i<nsnpbuf; i++) snp_id_map_buf.erase(adjsnps[i]);
    
    std::map<string,int>::iterator iter;
    for(iter=snp_id_map_buf.begin(); iter!=snp_id_map_buf.end(); iter++)
        snp_id_map.erase(iter->first);
    
    nsnpbuf = snp_id_map.size();
    remain_snp_indx.clear(); remain_snp_indx.resize(nsnpbuf);
    for(iter=snp_id_map.begin(), i=0; iter!=snp_id_map.end(); iter++, i++) remain_snp_indx[i] = iter->second;
    
    stable_sort(remain_snp_indx.begin(), remain_snp_indx.end());
}

void gcta::update_mtcojo_snp_rm(vector<string> adjsnps, map<string,int> &snp_id_map, vector<int> &remain_snp_indx) {
    
    int i=0, nsnpbuf=adjsnps.size();
    std::map<string,int>::iterator iter;
    
    for(i=0; i<nsnpbuf; i++) snp_id_map.erase(adjsnps[i]);
    
    nsnpbuf = snp_id_map.size();
    remain_snp_indx.clear(); remain_snp_indx.resize(nsnpbuf);
    for(iter=snp_id_map.begin(), i=0; iter!=snp_id_map.end(); iter++, i++) remain_snp_indx[i] = iter->second;

    stable_sort(remain_snp_indx.begin(), remain_snp_indx.end());
}

vector<string> gcta::read_snp_metafile_txt(string metafile, map<string,int> &gws_snp_name_map, double thresh) {
    ifstream meta_snp(metafile.c_str());
    if (!meta_snp)
         LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    
    string strbuf="", valbuf = "";
    vector<string> snplist;
    int line_number=0;

    // Read the summary data
    while(std::getline(meta_snp, strbuf)) {
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + metafile + "] should be in GCTA-COJO format, line " + to_string(line_number) + ".");
        }
        if(line_number==1) continue;
        // keep significant SNPs
        snplist.push_back(line_elements[0]);
        valbuf = line_elements[6];
        double pval_buf = 1.0;
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") pval_buf = atof(valbuf.c_str());
        if( pval_buf < thresh && (gws_snp_name_map.find(line_elements[0]) == gws_snp_name_map.end()) ) {
            int size = gws_snp_name_map.size();
            gws_snp_name_map.insert(pair<string,int>(line_elements[0], size));
        }
    }
    meta_snp.close();
    return snplist;
}

vector<string> gcta::read_snp_metafile_gz(string metafile, map<string,int> &gws_snp_name_map, double thresh) {
    string strbuf="", snpbuf="";
    vector<string> snplist;
    int line_number=0;
    const int MAX_LINE_LENGTH = 1024;
    char buf[MAX_LINE_LENGTH];
    
    gzifstream meta_snp(metafile.c_str());
    if (!meta_snp)
         LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    
    string err_msg = "Failed to read [" + metafile + "]. An error occurs in line ";

    // Read the summary data
    while(1) {
        meta_snp.getline(buf, MAX_LINE_LENGTH, '\n');
        if(meta_snp.fail() || !meta_snp.good()) break;
        line_number ++;
        if(line_number==1) continue;
        stringstream ss(buf);
        if (!(ss >> snpbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        // Save the SNPs
        snplist.push_back(snpbuf);
        int i = 0;
        for( i=0; i<6; i++)
            if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        // Save the p-value
        double pval_buf = 1.0;
        if(strbuf!="." && strbuf!="NA" && strbuf!="NAN") pval_buf = atof(strbuf.c_str());
        if( pval_buf < thresh && (gws_snp_name_map.find(snpbuf) == gws_snp_name_map.end()) ) {
            int size = gws_snp_name_map.size();
            gws_snp_name_map.insert(pair<string,int>(snpbuf, size));
        }
    }

    meta_snp.close();
    return snplist;
}

double gcta::read_single_metafile_txt(string metafile, map<string, int> id_map,
                         vector<string> &snp_a1, vector<string> &snp_a2,
                         eigenVector &snp_freq, eigenVector &snp_b,
                         eigenVector &snp_se, eigenVector &snp_pval,
                         eigenVector &snp_n, vector<bool> &snp_flag) {
   
    ifstream meta_raw(metafile.c_str());
    if (!meta_raw)
        LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    string strbuf="", valbuf="";
    int line_number=0, snp_indx=0;
    map<string, int>::iterator iter;
    // Read the summary data
    double pval_thresh = 0.5, h_buf = 0.0, vp_buf = 0.0, median_vp = 0.0;
    bool missing_flag = false;
    vector<double> vec_vp_buf;
     
    while(std::getline(meta_raw, strbuf)) {
        missing_flag = false;
        line_number++;
        std::istringstream linebuf(strbuf);
        std::istream_iterator<string> begin(linebuf), end;
        vector<string> line_elements(begin, end);
        if(line_elements.size() != 8) {
            LOGGER.e(0, "the GWAS summary data file [" + metafile + "] should be in GCTA-COJO format, line " + to_string(line_number) + ".");
        }        
        // Read the summary data
        if(line_number==1) continue;
        iter = id_map.find(line_elements[0]);      
        if(iter == id_map.end()) continue;     

        snp_indx = iter->second;         
        snp_a1[snp_indx] = line_elements[1]; snp_a2[snp_indx] = line_elements[2];
        StrFunc::to_upper(snp_a1[snp_indx]); StrFunc::to_upper(snp_a2[snp_indx]);
        valbuf =line_elements[3];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_freq(snp_indx) = atof(valbuf.c_str());
        else { snp_freq(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[4];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_b(snp_indx) = atof(valbuf.c_str());
        else { snp_b(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[5];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_se(snp_indx) = atof(valbuf.c_str());
        else { snp_se(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[6];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_pval(snp_indx) = atof(valbuf.c_str());
        else { snp_pval(snp_indx) = nan(""); missing_flag=true; }
        valbuf = line_elements[7];
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_n(snp_indx) = atof(valbuf.c_str());
        else { snp_n(snp_indx) = nan(""); missing_flag=true; }
        snp_flag[snp_indx] = true;

        if(!missing_flag) {
            h_buf = 2 * snp_freq(snp_indx) * ( 1- snp_freq(snp_indx) );
            vp_buf = h_buf * snp_b(snp_indx) * snp_b(snp_indx) + h_buf * snp_n(snp_indx) * snp_se(snp_indx) * snp_se(snp_indx);
            vec_vp_buf.push_back(vp_buf);
        }        
    }
    meta_raw.close();
    if(vec_vp_buf.size()>0) median_vp = CommFunc::median(vec_vp_buf);

    return median_vp;
}

double gcta::read_single_metafile_gz(string metafile, map<string, int> id_map,
                         vector<string> &snp_a1, vector<string> &snp_a2,
                         eigenVector &snp_freq, eigenVector &snp_b,
                         eigenVector &snp_se, eigenVector &snp_pval,
                         eigenVector &snp_n, vector<bool> &snp_flag) {
   
    const int MAX_LINE_LENGTH = 1024;
    char buf[MAX_LINE_LENGTH];

    gzifstream meta_raw(metafile.c_str());
    if (!meta_raw)
        LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");
    string err_msg = "Failed to read [" + metafile + "]. An error occurs in line ";

    int line_number=0, snp_indx=0;
    map<string, int>::iterator iter;
    // Read the summary data
    double pval_thresh = 0.5, h_buf = 0.0, vp_buf = 0.0, median_vp = 0.0;
    bool missing_flag = false;
    string strbuf = "";
    vector<double> vec_vp_buf;
    while(1) {
        meta_raw.getline(buf, MAX_LINE_LENGTH, '\n');
        if(meta_raw.fail() || !meta_raw.good()) break;
        missing_flag = false;
        line_number++;
        if(line_number==1) continue;
         
        string snpbuf = "", strbuf = "", valbuf = "";
        stringstream ss(buf);
        // SNP
        if (!(ss >> snpbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        iter = id_map.find(snpbuf);
        if(iter == id_map.end()) continue;     
        snp_indx = iter->second;
        // a1 
        if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        StrFunc::to_upper(strbuf);
        snp_a1[snp_indx] = strbuf;
        // a2
        if (!(ss >> strbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        StrFunc::to_upper(strbuf);
        snp_a2[snp_indx] = strbuf;
        // freq
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_freq(snp_indx) = atof(valbuf.c_str());
        else { snp_freq(snp_indx) = nan(""); missing_flag=true; }
        // b
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_b(snp_indx) = atof(valbuf.c_str());
        else { snp_b(snp_indx) = nan(""); missing_flag=true; }
        // se
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_se(snp_indx) = atof(valbuf.c_str());
        else { snp_se(snp_indx) = nan(""); missing_flag=true; }
        // p-value
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_pval(snp_indx) = atof(valbuf.c_str());
        else { snp_pval(snp_indx) = nan(""); missing_flag=true; }
        // n
        if (!(ss >> valbuf)) LOGGER.e(0, err_msg + to_string(line_number) + "," + buf);
        if(valbuf!="." && valbuf!="NA" && valbuf!="NAN") snp_n(snp_indx) = atof(valbuf.c_str());
        else { snp_n(snp_indx) = nan(""); missing_flag=true; }
        snp_flag[snp_indx] = true;

        if(!missing_flag) {
            h_buf = 2 * snp_freq(snp_indx) * ( 1- snp_freq(snp_indx) );
            vp_buf = h_buf * snp_b(snp_indx) * snp_b(snp_indx) + h_buf * snp_n(snp_indx) * snp_se(snp_indx) * snp_se(snp_indx);
            vec_vp_buf.push_back(vp_buf);
        }        
    }
    meta_raw.close();
    if(vec_vp_buf.size()>0) median_vp = CommFunc::median(vec_vp_buf);

    return median_vp;
}

vector<string> gcta::remove_bad_snps(vector<string> snp_name, vector<int> snp_remain, vector<vector<bool>> snp_flag, vector<vector<string>> &snp_a1, vector<vector<string>> &snp_a2, eigenMatrix &snp_freq,  
                                    eigenMatrix &snp_b, eigenMatrix snp_se, eigenMatrix snp_pval, eigenMatrix snp_n, map<string,int> plink_snp_name_map, vector<string> snp_ref_a1, vector<string> snp_ref_a2, 
                                    int ntarget, int ncovar, string outfile_name) {
    int i=0, j=0, nsnp = snp_remain.size(), npheno = ntarget+ncovar;
    double snp_freq_bak = 0.0;
    string snp_a1_bak = "", snp_a2_bak = ""; 
    vector<string> badsnps;
    vector<int> bad_indx;
    map<string,int>::iterator iter;
    vector<vector<bool>> flip_flag(npheno);

    for(i=0; i<npheno; i++) {
        flip_flag[i].resize(nsnp);
        for(j=0; j<nsnp; j++) flip_flag[i][j] = false;
    }

    for( i=0; i<nsnp; i++) {
        bool iterFlag=true;
        int snpindx = -9;
        vector<string> allelebuf;
        // Choose a SNP
        iter = plink_snp_name_map.find(snp_name[snp_remain[i]]);
        if(iter!=plink_snp_name_map.end()) snpindx = iter->second;
        // Alleles
        if(snpindx > 0) {
            allelebuf.push_back(snp_ref_a1[snpindx]); allelebuf.push_back(snp_ref_a2[snpindx]);
        } else {
            for(j=0; j<npheno; j++) {
                if(snp_flag[j][snp_remain[i]]) {
                    allelebuf.push_back(snp_a1[j][snp_remain[i]]); 
                    allelebuf.push_back(snp_a2[j][snp_remain[i]]);
                    break;
                }
            }
        }

        snp_a1_bak = allelebuf[0]; snp_a2_bak = allelebuf[1];
        for(j=0; j<npheno; j++) {
            if(!snp_flag[j][snp_remain[i]]) continue;
            // Alleles
            allelebuf.push_back(snp_a1[j][snp_remain[i]]);
            allelebuf.push_back(snp_a2[j][snp_remain[i]]);
            // Removing SNPs with missing value
            if( std::isnan(snp_b(snp_remain[i],j)) || std::isnan(snp_se(snp_remain[i],j)) || std::isnan(snp_pval(snp_remain[i],j)) || std::isnan(snp_n(snp_remain[i],j)) ) {
                iterFlag=false; break;
            }
            // Removing SNPs with extremely small SE
            if(snp_se(snp_remain[i],j) < eps) {
                iterFlag=false; break;
            }
            // Alignment of b
            if(iterFlag && (allelebuf[0] != snp_a1[j][snp_remain[i]])) {
                snp_b(snp_remain[i],j) = -1*snp_b(snp_remain[i],j);
                if(!std::isnan(snp_freq(snp_remain[i],j))) 
                    snp_freq(snp_remain[i],j)  =  1 - snp_freq(snp_remain[i],j);
                std::swap(snp_a1[j][snp_remain[i]], snp_a2[j][snp_remain[i]]);
            }
            snp_freq_bak = snp_freq(snp_remain[i],j);
        }

        // Remove SNPs with multiple alleles
        stable_sort(allelebuf.begin(), allelebuf.end());      
        allelebuf.erase(unique(allelebuf.begin(), allelebuf.end()), allelebuf.end());
        if(allelebuf.size()!=2) iterFlag=false;
        // Collect bad SNPs
        if(iterFlag==false) {
            badsnps.push_back(snp_name[snp_remain[i]]);
        }
        
        // update the estimates for the first phenotype
        if(!snp_flag[0][snp_remain[i]]) {
            snp_a1[0][snp_remain[i]] = snp_a1_bak;
            snp_a2[0][snp_remain[i]] = snp_a2_bak;
            snp_freq(snp_remain[i],0) = snp_freq_bak;
        }
    }

    if (!badsnps.empty()) {
        string badsnpfile = outfile_name + ".badsnps", strbuf="";
        ofstream obadsnp(badsnpfile.c_str());
        if(!obadsnp) LOGGER.e(0, "cannot open file [" + badsnpfile + "] to write bad SNPs.");
        int nbadsnps = badsnps.size();
        for (i = 0; i < nbadsnps; i++) 
            obadsnp << badsnps[i] << endl;
        obadsnp.close();
        LOGGER.i(0,  to_string(nbadsnps) + " SNPs have missing value or mismatched alleles. These SNPs have been saved in [" + badsnpfile + "].");
    }
   
    stable_sort(badsnps.begin(), badsnps.end());
    return badsnps;
}

vector<string> gcta::remove_freq_diff_snps(vector<string> meta_snp_name, vector<int> meta_snp_remain, map<string,int> snp_name_map, vector<double> ref_freq, eigenMatrix meta_freq, vector<vector<bool>> snp_flag, int ntrait, double freq_thresh, string outfile_name) {
    int i = 0, nsnp = meta_snp_remain.size(), nsnp_ttl = meta_snp_remain.size();
    string snpbuf="";
    vector<string> afsnps;
    map<string,int>::iterator iter_ref;

    for( i=0; i<nsnp; i++ ) {
        int refsnp_index = 0;
        vector<double> snpfreq;

        // AF from the reference sample
        snpbuf = meta_snp_name[meta_snp_remain[i]];
        iter_ref = snp_name_map.find(snpbuf);
        if( iter_ref != snp_name_map.end()) {
            refsnp_index = iter_ref -> second;
            snpfreq.push_back(ref_freq[refsnp_index]/2);
        }

        // AF from the GWAS summary
        int j = 0;
        for( j=0; j<ntrait; j++ ) {
            if(!snp_flag[j][meta_snp_remain[i]]) continue;
            if(std::isnan(meta_freq(meta_snp_remain[i],j))) continue;
            snpfreq.push_back(meta_freq(meta_snp_remain[i],j));
        }

        // Test the difference
        int t1 = 0, t2 = 0, nsnp_test = snpfreq.size();
        double max_freq_dev = 0.0, freq_dev = 0.0;
        for( t1=0; t1<nsnp_test-1; t1++) {
            for(t2=t1+1; t2<nsnp_test; t2++) {
                freq_dev = fabs(snpfreq[t1]-snpfreq[t2]);
                if(max_freq_dev < freq_dev) max_freq_dev = freq_dev;
            }
        }
        if(max_freq_dev > freq_thresh) afsnps.push_back(snpbuf);     
    }  

    if (!afsnps.empty()) {
        string afsnpfile = outfile_name + ".freq.badsnps", strbuf="";
        ofstream oafsnp(afsnpfile.c_str());
        if(!oafsnp) LOGGER.e(0, "cannot open file [" + afsnpfile + "] to write bad SNPs.");
        int nafsnps = afsnps.size();
        for (i = 0; i < nafsnps; i++) oafsnp << afsnps[i] << endl;
        oafsnp.close();
        LOGGER.i(0,  to_string(nafsnps) + " SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample. These SNPs have been saved in [" + afsnpfile + "].");
        if(nafsnps > nsnp_ttl*0.05) 
            LOGGER.e(0, "there are too many SNPs that have large difference in allele frequency. Please check the GWAS summary data.");
    }
    return(afsnps);
}

vector<string> gcta::remove_mono_snps(map<string,int> snp_name_map, vector<double> ref_snpfreq, string outfile_name) {
    int i = 0, n_raresnp = 0;
    vector<string> afsnps;
    map<string,int>::iterator iter;

    for(iter=snp_name_map.begin(); iter!=snp_name_map.end(); iter++) {
        double af = ref_snpfreq[iter->second]/2;
        if(af > 0.5) af = 1.0 - af;
        if(CommFunc::FloatEqual(af, 0.0)) afsnps.push_back(iter->first);
        if(af < 0.01) n_raresnp++;
    }

    if(n_raresnp > 0) LOGGER.w(0,  "There are " + to_string(n_raresnp) + " SNPs with MAF < 0.01 in the reference sample.");

    if (!afsnps.empty()) {
        string afsnpfile = outfile_name + ".mono.badsnps", strbuf="";
        ofstream oafsnp(afsnpfile.c_str());
        if(!oafsnp) LOGGER.e(0, "cannot open file [" + afsnpfile + "] to write bad SNPs.");
        int nafsnps = afsnps.size();
        for (i = 0; i < nafsnps; i++) oafsnp << afsnps[i] << endl;
        oafsnp.close();
        LOGGER.i(0,  to_string(nafsnps) + " monomorphic SNP(s) have been saved in [" + afsnpfile + "].");
    }

    return(afsnps);
}

vector<string> gcta::filter_meta_snp_pval(vector<string> snp_name, vector<int> remain_snp_indx,  eigenMatrix snp_pval, int start_indx, int end_indx, vector<vector<bool>> snp_flag, double pval_thresh) {
    
    int i=0, j=0, nsnpbuf = remain_snp_indx.size();
    vector<string> kpsnps;
   
    for(i=0; i<nsnpbuf; i++) {
        for(j=start_indx; j<end_indx; j++) {
            if((snp_pval(remain_snp_indx[i],j) <  pval_thresh) && (snp_flag[j][remain_snp_indx[i]])) {
                kpsnps.push_back(snp_name[remain_snp_indx[i]]);
                break;
            }
        }
    }
    return kpsnps;
}
