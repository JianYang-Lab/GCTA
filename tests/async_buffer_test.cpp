#include "infra/AsyncBuffer.hpp"

#include <cstdint>
#include <iostream>
#include <thread>

int main(){
    const uint64_t block_count = 200;
    const uint64_t block_size = 8;
    for(int round = 0; round < 3; round++){
        AsyncBuffer<uint64_t> buffer(block_size);
        if(!buffer.init_status()) return 1;
        std::thread producer([&](){
            for(uint64_t block = 0; block < block_count; block++){
                uint64_t *data = buffer.start_write();
                for(uint64_t index = 0; index < block_size; index++){
                    data[index] = block * block_size + index;
                }
                if(block + 1 == block_count) buffer.setEOF();
                buffer.end_write();
                if(block % 17 == 0) std::this_thread::yield();
            }
        });
        for(uint64_t block = 0; block < block_count; block++){
            auto [data, eof] = buffer.start_read();
            for(uint64_t index = 0; index < block_size; index++){
                if(data[index] != block * block_size + index){
                    std::cerr << "data mismatch at block " << block << " index " << index << std::endl;
                    producer.join();
                    return 1;
                }
            }
            if(eof != (block + 1 == block_count)){
                std::cerr << "EOF mismatch at block " << block << std::endl;
                producer.join();
                return 1;
            }
            buffer.end_read();
            if(block % 19 == 0) std::this_thread::yield();
        }
        producer.join();
    }
    return 0;
}
