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
    Timer total_seed_selection;
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

    std::vector<TimerAggregator> subTimers;

    TimerAggregator(){}
    ~TimerAggregator(){}

    nlohmann::json buildPathTimesJson(const int world_size) {
        nlohmann::json executionPathTimes;

        for (const auto & e : subTimers) {
            nlohmann::json subTime{
                {"MaxKCover", this->getWorldTimes(world_size, e.max_k_globalTimer.resolveTimer() + e.max_k_localTimer.resolveTimer())},  
                {"MPI_Gather", this->getWorldTimes(world_size, e.allGatherTimer.resolveTimer())}, 
            };

            executionPathTimes.push_back(subTime);
        }

        return executionPathTimes;
    }

    nlohmann::json buildLazyLazyTimeJson(const int world_size, const double total_runtime) {
        nlohmann::json timeReport{
            {"Sampling", this->getWorldTimes(world_size, this->samplingTimer.resolveTimer())}, 
            {"MPI_AllToAll", this->getWorldTimes(world_size, this->allToAllTimer.resolveTimer())}, 
            {"GranularExecutionPaths", this->buildPathTimesJson(world_size)},
            {"MPI_Broadcast", this->getWorldTimes(world_size, this->broadcastTimer.resolveTimer())},
            {"TotalSeedSelection", this->getWorldTimes(world_size, this->total_seed_selection.resolveTimer())},
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
            {"TotalSeedSelection", this->getWorldTimes(world_size, this->total_seed_selection.resolveTimer())},
            {"Total", this->getTotalRuntimes(world_size, total_runtime)} 
        };

        return timeReport;
    }
};

