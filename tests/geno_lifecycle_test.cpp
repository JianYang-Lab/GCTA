#include "data/Geno.h"
#include "data/Marker.h"
#include "data/Pheno.h"

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#ifdef __linux__
size_t open_file_count(){
    return std::distance(
        std::filesystem::directory_iterator("/proc/self/fd"),
        std::filesystem::directory_iterator());
}
#else
size_t open_file_count(){
    return 0;
}
#endif

bool check_layout(uint32_t sample_count){
    constexpr uint32_t marker_count = 2;
    const uint32_t bytes_per_marker = (sample_count + 3) / 4;
    const uint32_t words_per_marker = (bytes_per_marker + 7) / 8;
    std::vector<uint8_t> input(marker_count * bytes_per_marker + 8 * sizeof(uint64_t));
    for(size_t index = 0; index < input.size(); index++){
        input[index] = static_cast<uint8_t>(11 + 37 * index);
    }
    std::vector<uint64_t> keep((sample_count + 63) / 64, 0);
    for(uint32_t index = 0; index < sample_count; index++){
        keep[index / 64] |= 1ULL << (index % 64);
    }
    std::vector<uint64_t> actual(marker_count * words_per_marker, 0);
    std::vector<uint64_t> expected(marker_count * words_per_marker, 0);
    for(uint32_t marker = 0; marker < marker_count; marker++){
        memcpy(
            reinterpret_cast<uint8_t *>(expected.data()) + marker * words_per_marker * sizeof(uint64_t),
            input.data() + marker * bytes_per_marker,
            bytes_per_marker);
    }
    Geno::move_geno(
        input.data(),
        keep.data(),
        sample_count,
        sample_count,
        marker_count,
        actual.data());
    if(actual == expected) return true;
    std::cerr << "stale genotype layout for " << sample_count << " samples" << std::endl;
    return false;
}

int main(int argc, char **argv){
    if(argc != 2) return 2;
    if(!check_layout(8) || !check_layout(16) || !check_layout(40)) return 1;

    std::map<std::string, std::vector<std::string>> options;
    options["--bfile"] = {argv[1]};
    options["--out"] = {"geno_lifecycle"};
    options["--make-bed"] = {};
    Pheno::registerOption(options);
    Marker::registerOption(options);
    Geno::registerOption(options);

    const size_t initial_files = open_file_count();
    for(int round = 0; round < 4; round++){
        std::vector<uint64_t> packed;
        size_t genotype_bytes = 0;
        {
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.init_keep();
            packed.assign((pheno.count_keep() + 31) / 32, 0x0123456789abcdefULL + round);
            genotype_bytes = (pheno.count_keep() + 3) / 4;
            geno.save_bed(packed.data(), 1);
            geno.closeOut();
        }
        std::ifstream stream("geno_lifecycle.bed", std::ios::binary);
        const std::string written((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
        const std::string header("\x6c\x1b\x01", 3);
        if(written.size() != genotype_bytes + 3 || written.compare(0, 3, header) != 0 ||
            memcmp(written.data() + 3, packed.data(), genotype_bytes) != 0){
            std::cerr << "repeated BED output reused stale file state" << std::endl;
            return 1;
        }
    }
    const size_t final_files = open_file_count();
    if(final_files != initial_files){
        std::cerr << "genotype file handles leaked: " << initial_files
                  << " before, " << final_files << " after" << std::endl;
        return 1;
    }
    return 0;
}
