
typedef struct {
    std::map<int, std::vector<int>> allKMSeeds,
    std::vector<std::pair<unsigned int, std::vector<unsigned int>>> localCandidateSets
} SolutionCandidateSets ;

class ApproximatorGroup {
    protected: 
    const MPI_Comm groupWorld;
    const std::vector<int> &vertexToProcess;

    ApproximatorGroup(
        const std::vector<int> &vertexToProcess
    ) : vertexToProcess(vertexToProcess) {

    }

    public:
    virtual SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &aggregateSets,
        std::pair<std::vector<unsigned int>, ssize_t> &localCandidateSet,
        const int kprime, 
        const size_t theta
    ) const = 0;
}

class DummyApproximatorGroup : public ApproximatorGroup {
    public:
    DummyApproximatorStrategy(const MPI_Comm groupWorld) : groupWorld(groupWorld) {}

    SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &aggregateSets,
        std::pair<std::vector<unsigned int>, ssize_t> &localCandidateSet,
        const int kprime, 
        const size_t theta
    ) const override {
        int group_rank; 
        MPI_Comm_rank(groupWorld, group_rank);
        
        if (group_rank == 0) {
            return {
                std::map<int, std::vector<int>>{
                    0, std::vector<unsigned int>{0,1,2,3,4,5,6,7,8,9},
                    1, std::vector<unsigned int>{0,1,2,3,4,5,6,7,8,9},
                    2, std::vector<unsigned int>{0,1,2,3,4,5,6,7,8,9}
                },
                std::vector<std::pair<unsigned int, std::vector<unsigned int>>>{
                    80000, std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
                }
            };
        } else {
            return  {
                std::map<int, std::vector<int>>(),
                std::vector<std::pair<unsigned int, std::vector<unsigned int>>>()
            }
        }
    }
};

template <typename GraphTy, typename ConfTy>
class LazyLazyApproximatorGroup : public ApproximatorGroup {
    private:
    TimerAggregator &timeAggregator;
    const ConfTy &CFG;
    const CommunicationEngine<GraphTy> &cEngine;

    static std::pair<std::vector<unsigned int>, ssize_t> SolveKCover(
        const int k,
        const int kprime,
        const size_t max_element, 
        TimerAggregator &timeAggregator,
        const std::map<int, std::vector<int>>& elements
    )
    {
        MaxKCover<GraphTy> localKCoverEngine(k, kprime, max_element, timeAggregator);
        localKCoverEngine.useLazyGreedy();
        return localKCoverEngine.run_max_k_cover(elements);
    }

    public:
    LazyLazyApproximatorGroup(
        TimerAggregator &timeAggregator,
        const ConfTy &input_CFG,
        const CommunicationEngine<GraphTy> &cEngine,
    ) : timeAggregator(timeAggregator), CFG(input_CFG), cEngine(cEngine) {}

    SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &aggregateSets,
        std::pair<std::vector<unsigned int>, ssize_t> &localCandidateSet,
        const int kprime, 
        const size_t theta
    ) const override {
        if (localCandidateSet.first.size() == 0) {
            this->timeAggregator.max_k_localTimer.startTimer();

            std::pair<std::vector<unsigned int>, ssize_t> newSeeds = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, aggregateSets
            );

            localCandidateSet.first = newSeeds.first;
            localCandidateSet.second = newSeeds.second;

            this->timeAggregator.max_k_localTimer.endTimer();
        }

        // spdlog::get("console")->info("all gather...");
        this->timeAggregator.allGatherTimer.startTimer();

        std::vector<unsigned int> globalAggregation;
        size_t totalData = this->cEngine.GatherPartialSolutions(localCandidateSet, aggregateSets, globalAggregation);

        this->timeAggregator.allGatherTimer.endTimer();

        std::pair<std::vector<unsigned int>, size_t> approximated_solution;
        std::map<int, std::vector<int>> bestKMSeeds;
        std::vector<std::pair<unsigned int, std::vector<unsigned int>>> candidateSets;

        if (this->cEngine.GetRank() == 0) {

            // spdlog::get("console")->info("unpacking local seeds in global process...");

            this->timeAggregator.allGatherTimer.startTimer();
            candidateSets = this->cEngine.aggregateLocalKSeeds(bestKMSeeds, globalAggregation.data(), totalData);

            this->timeAggregator.allGatherTimer.endTimer();

            // spdlog::get("console")->info("global max_k_cover...");
            this->timeAggregator.max_k_globalTimer.startTimer();

            approximated_solution = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, bestKMSeeds
            );

            candidateSets.push_back(std::make_pair(
                static_cast<unsigned int>(approximated_solution.second),
                approximated_solution.first
            ));

            this->timeAggregator.max_k_globalTimer.endTimer();
        }

        return {
            bestKMSeeds,
            candidateSets
        };
    }
};