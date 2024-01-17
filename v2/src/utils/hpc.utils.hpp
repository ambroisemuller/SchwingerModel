


namespace HPC {
    
    void MultithreadedLoop(std::function<void(int, int, int)> func, int start, int end) {
        int num_threads = std::thread::hardware_concurrency();
        
        std::vector<std::thread> threads;
        int range = (end - start) / num_threads;

        for (int t = 0; t < num_threads; ++t) {
            int start_idx = start + t * range;
            int end_idx = (t == num_threads - 1) ? end : start_idx + range;

            threads.emplace_back([=, &func]() {
                func(start_idx, end_idx, t);
            });
        }

        for (auto &thread : threads) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
}
