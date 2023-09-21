#include <vector>

#include <utility>
#include <chrono>

class Timer 
{
    private:
    std::vector<std::chrono::duration<double, std::milli>> timer;
    std::chrono::_V2::system_clock::time_point start;
    std::chrono::_V2::system_clock::time_point end;

    public:
    Timer(){}
    ~Timer(){}

    void startTimer()
    {
        start = std::chrono::high_resolution_clock::now();
    }

    void endTimer()
    {
        end = std::chrono::high_resolution_clock::now();
        timer.push_back(end - start);
    }

    double resolveTimerDEBUG() const
    {
        double totalTime = 0;
        for (const auto time : timer)
        {
            std::cout << time.count() << std::endl;
            totalTime += time.count();
        }
        return totalTime;
    }

    double resolveTimer() const
    {
        double totalTime = 0;
        for (const auto time : timer)
        {
            totalTime += time.count();
        }
        return totalTime;
    }
};

class TimerAggregator
{
    private:
    static std::vector<double> getWorldTimes(const int worldSize, const double localTime) {
        std::vector<double> res(worldSize, 0);
        MPI_Gather(&localTime, 1, MPI_DOUBLE, res.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return res;
    }

    static std::vector<double> getTotalRuntimes(const int worldSize, const double totalRuntime) {
        std::vector<double> res(worldSize, 0);
        MPI_Gather(&totalRuntime, 1, MPI_DOUBLE, res.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return res;
    }

    public:
    Timer samplingTimer;
    Timer max_k_localTimer;
    Timer max_k_globalTimer;
    Timer allGatherTimer;
    Timer allToAllTimer;
    Timer broadcastTimer;
    Timer totalGlobalStreamTimer;
    Timer localStreamTimer;
    Timer sendTimer;
    Timer receiveTimer;
    Timer initBucketTimer;
    Timer insertIntoBucketTimer;

    Timer processingReceiveTimer;
    Timer atomicUpdateTimer;

    Timer totalSendTimer;
    Timer totalReceiveTimer;

    TimerAggregator(){}
    ~TimerAggregator(){}

    nlohmann::json buildLazyLazyTimeJson(const int world_size, const double total_runtime) {
        nlohmann::json timeReport{
            {"Sampling", this->getWorldTimes(world_size, this->samplingTimer.resolveTimer())}, 
            {"MPI_AllToAll", this->getWorldTimes(world_size, this->allToAllTimer.resolveTimer())}, 
            {"MaxKCover", this->getWorldTimes(world_size, this->max_k_globalTimer.resolveTimer() + this->max_k_localTimer.resolveTimer())},  
            {"MPI_Gather", this->getWorldTimes(world_size, this->allGatherTimer.resolveTimer())}, 
            {"MPI_Broadcast", this->getWorldTimes(world_size, this->broadcastTimer.resolveTimer())},
            {"Total", this->getTotalRuntimes(world_size, total_runtime)} 
        };

        return timeReport;
    }

    nlohmann::json buildStreamingTimeJson(const int world_size, const double total_runtime) {
        nlohmann::json timeReport{
            {"Sampling", this->getWorldTimes(world_size, this->samplingTimer.resolveTimer())}, 
            {"MPI_AllToAll", this->getWorldTimes(world_size, this->allToAllTimer.resolveTimer())}, 
            {"SelectNextSeedUsingMaxKCover", this->getWorldTimes(world_size, this->max_k_localTimer.resolveTimer())},  
            {"SendingNextSeed", this->getWorldTimes(world_size, this->sendTimer.resolveTimer())},
            {"TotalSendTime", this->getWorldTimes(world_size, this->totalSendTimer.resolveTimer())},
            {"InitializingBuckets_ListeningThread", this->getWorldTimes(world_size, this->initBucketTimer.resolveTimer())},
            {"WaitingToReceiveNextSeed_ListeningThread", this->getWorldTimes(world_size, this->receiveTimer.resolveTimer())},
            {"InsertingSeed_BucketingThreads", this->getWorldTimes(world_size, this->max_k_globalTimer.resolveTimer())},
            {"ProcessingSentSeedForBuckets_ListeningThread", this->getWorldTimes(world_size, this->processingReceiveTimer.resolveTimer())},
            {"AtomicUpdateForBuckets_ListeningThread", this->atomicUpdateTimer.resolveTimer()},
            {"TotalGlobalStreamingTime", this->totalGlobalStreamTimer.resolveTimer()},
            {"MPI_Broadcast", this->getWorldTimes(world_size, this->broadcastTimer.resolveTimer())},
            {"Total", this->getTotalRuntimes(world_size, total_runtime)} 
        };

        return timeReport;
    }


    // void OutputDiagnosticData()
    // {
    //     if (CFG.use_streaming == true)
    //     {
    //         std::cout << " --- SHARED --- " << std::endl; 
    //         std::cout << "Samping time: " << this->samplingTimer.resolveTimer() << std::endl;
    //         std::cout << "AlltoAll time: " << this->allToAllTimer.resolveTimer() << std::endl;
    //         std::cout << "Receive Broadcast: " << this->broadcastTimer.resolveTimer() << std::endl;

    //         std::cout << " --- SENDER --- " << std::endl; 
    //         std::cout << "Select Next Seed: " << this->max_k_localTimer.resolveTimer() << std::endl;
    //         std::cout << "Send Next Seed: " << this->sendTimer.resolveTimer() << std::endl;
    //         std::cout << "Total Send Time: " << this->totalSendTimer.resolveTimer() << std::endl;
            
    //         std::cout << " --- RECEIVER --- " << std::endl; 
    //         std::cout << "Initialize Buckets: " << this->initBucketTimer.resolveTimer() << std::endl;
    //         std::cout << "Receive Next Seed: " << this->receiveTimer.resolveTimer() << std::endl;
    //         std::cout << "Insert Into Buckets: " << this->max_k_globalTimer.resolveTimer() << std::endl;
    //         std::cout << "Handling received data (inserting into matrix and copying from buffer): " << this->processingReceiveTimer.resolveTimer() << std::endl; 
    //         std::cout << "Atomic Update (receiver side): " << this->atomicUpdateTimer.resolveTimer() << std::endl; 
    //         std::cout << "Total Global Streaming Time: " << this->totalGlobalStreamTimer.resolveTimer() << std::endl;
    //     } else {
    //         std::cout << " --- SHARED --- " << std::endl; 
    //         std::cout << "Samping time: " << this->samplingTimer.resolveTimer() << std::endl;
    //         std::cout << "f score Broadcast time: " << this->broadcastTimer.resolveTimer() << std::endl;
    //         std::cout << "AlltoAll time: " << this->allToAllTimer.resolveTimer() << std::endl;
    //         std::cout << "AllGather time: " << this->allGatherTimer.resolveTimer() << std::endl;

    //         std::cout << " --- LOCAL --- " << std::endl; 
    //         std::cout << "Local max-cover time: " << this->max_k_localTimer.resolveTimer() << std::endl;

    //         std::cout << " --- GLOBAL --- " << std::endl; 
    //         std::cout << "Global max-cover time: " << this->max_k_globalTimer.resolveTimer() << std::endl;
    //     }
    // }

};

