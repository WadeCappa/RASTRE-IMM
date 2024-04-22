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
    Timer totalSeedSelectionTimer;
    Timer sendTimer;
    Timer receiveTimer;
    Timer initBucketTimer;
    Timer insertIntoBucketTimer;

    Timer processingReceiveTimer;
    Timer atomicUpdateTimer;

    Timer totalReceiveTimer;

    Timer broadcastSeeds_OPIMC;
    
    Timer countUncoveredR1_OPIMC;
    Timer reduceUncoveredR1_OPIMC;
    Timer selectNextBestKR1_OPIMC;

    Timer countCoveredR2_OPIMC;
    Timer reduceCoveredR2_OPIMC;

    Timer barrierTime;

    TimerAggregator(){}
    ~TimerAggregator(){}

    nlohmann::json buildOpimcTimeJson(const int world_size) {
        nlohmann::json timeReport{
            {"Broadcasting Seeds", this->getWorldTimes(world_size, this->broadcastSeeds_OPIMC.resolveTimer())},
            {"Estimating alpha upper bound", nlohmann::json {
                {"Count uncovered on R1", this->getWorldTimes(world_size, this->countUncoveredR1_OPIMC.resolveTimer())},
                {"Reduce unconvered on R1", this->getWorldTimes(world_size, this->reduceUncoveredR1_OPIMC.resolveTimer())},
                {"Select next best k on R1", this->getWorldTimes(world_size, this->selectNextBestKR1_OPIMC.resolveTimer())}
            }},
            {"Getting influence on R2", nlohmann::json {
                {"Count covered on R2", this->getWorldTimes(world_size, this->countCoveredR2_OPIMC.resolveTimer())},
                {"Reduce covered on R2", this->getWorldTimes(world_size, this->reduceCoveredR2_OPIMC.resolveTimer())}
            }}
        };

        return timeReport;
    }

    nlohmann::json buildLazyLazyTimeJson(const int world_size, const double total_runtime) {
        nlohmann::json timeReport{
            {"Sampling", this->getWorldTimes(world_size, this->samplingTimer.resolveTimer())}, 
            {"MPI_AllToAll", this->getWorldTimes(world_size, this->allToAllTimer.resolveTimer())}, 
            {"Local MaxKCover", this->getWorldTimes(world_size, this->max_k_localTimer.resolveTimer())},  
            {"Global MaxKCover", this->getWorldTimes(world_size, this->max_k_globalTimer.resolveTimer())},  
            {"MPI_Gather", this->getWorldTimes(world_size, this->allGatherTimer.resolveTimer())}, 
            {"MPI_Broadcast", this->getWorldTimes(world_size, this->broadcastTimer.resolveTimer())},
            {"TotalSeedSelectionTime", this->getWorldTimes(world_size, this->totalSeedSelectionTimer.resolveTimer())},
            {"Time spent waiting after decoding AllToAll", this->getWorldTimes(world_size, this->barrierTime.resolveTimer())},
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
            {"InitializingBuckets_ListeningThread", this->initBucketTimer.resolveTimer()},
            {"WaitingToReceiveNextSeed_ListeningThread", this->receiveTimer.resolveTimer()},
            {"InsertingSeed_BucketingThreads", this->max_k_globalTimer.resolveTimer()},
            {"ProcessingSentSeedForBuckets_ListeningThread", this->processingReceiveTimer.resolveTimer()},
            {"AtomicUpdateForBuckets_ListeningThread", this->atomicUpdateTimer.resolveTimer()},
            {"TotalSeedSelectionTime", this->getWorldTimes(world_size, this->totalSeedSelectionTimer.resolveTimer())},
            {"MPI_Broadcast", this->getWorldTimes(world_size, this->broadcastTimer.resolveTimer())},
            {"Total", this->getTotalRuntimes(world_size, total_runtime)} 
        };

        return timeReport;
    }
};

