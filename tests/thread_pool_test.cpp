#include "infrastructure/ThreadPool.h"
#include <atomic>
#include <iostream>
#include <thread>

int main(){
    ThreadPool *pool = ThreadPool::GetPool(4);
    if(pool != ThreadPool::GetPool() || pool->getThreadCount() != 4) return 1;

    constexpr int job_count = 64;
    std::atomic<bool> start{false};
    std::atomic<int> completed{0};
    for(int index = 0; index < job_count; index++){
        pool->AddJob([&](){
            while(!start.load(std::memory_order_acquire)) std::this_thread::yield();
            completed.fetch_add(1, std::memory_order_relaxed);
        });
    }
    start.store(true, std::memory_order_release);
    pool->WaitAll();
    if(completed.load(std::memory_order_relaxed) != job_count || pool->GetRemainCount() != 0){
        std::cerr << "thread pool did not complete every job" << std::endl;
        return 1;
    }
    return 0;
}
