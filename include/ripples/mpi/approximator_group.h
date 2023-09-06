
typedef struct solutionCandidateSets {
    std::map<int, std::vector<int>> allKMSeeds;
    std::vector<std::pair<unsigned int, std::vector<unsigned int>>> localCandidateSets;
} SolutionCandidateSets ;

class ApproximatorGroup {
    protected: 
    const MPI_Comm groupWorld;
    const std::vector<int> &vertexToProcess;

    ApproximatorGroup(
        const MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess
    ) : vertexToProcess(vertexToProcess), groupWorld(groupWorld) {

    }

    public:
    virtual SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &localSolutionSpace,
        const std::pair<std::vector<unsigned int>, unsigned int> &previousGroupSolution,
        const int kprime, 
        const size_t theta
    ) const {
        std::cout << "entered the approximate function of the ApproximatorGroup base class, killing process" << std::endl;
        exit(1);
    };
};

class DummyApproximatorGroup : public ApproximatorGroup {
    public:
    DummyApproximatorGroup(
        const MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess
    ) : ApproximatorGroup(groupWorld, vertexToProcess) {}

    SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &localSolutionSpace,
        const std::pair<std::vector<unsigned int>, unsigned int> &previousGroupSolution,
        const int kprime, 
        const size_t theta
    ) const override {
        int group_rank; 
        MPI_Comm_rank(groupWorld, &group_rank);
        
        if (group_rank == 0) {
            return {
                std::map<int, std::vector<int>>{
                    {0, std::vector<int>{0,1,2,3,4,5,6,7,8,9}},
                    {1, std::vector<int>{0,1,2,3,4,5,6,7,8,9}},
                    {2, std::vector<int>{0,1,2,3,4,5,6,7,8,9}}
                },
                std::vector<std::pair<unsigned int, std::vector<unsigned int>>>{
                    std::make_pair(80000, std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9})
                }
            };
        } else {
            return  {
                std::map<int, std::vector<int>>(),
                std::vector<std::pair<unsigned int, std::vector<unsigned int>>>()
            };
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
        MPI_Comm group_world,
        const std::vector<int> &vertexToProcess,
        TimerAggregator &timeAggregator,
        const ConfTy &input_CFG,
        const CommunicationEngine<GraphTy> &cEngine
    ) : timeAggregator(timeAggregator), CFG(input_CFG), cEngine(cEngine), ApproximatorGroup(groupWorld, vertexToProcess) {}

    SolutionCandidateSets approximate(
        const std::map<int, std::vector<int>> &localSolutionSpace,
        const std::pair<std::vector<unsigned int>, unsigned int> &previousGroupSolution,
        const int kprime, 
        const size_t theta
    ) const override {
        auto workingCandidateSet = previousGroupSolution;
        if (workingCandidateSet.first.size() == 0) {
            this->timeAggregator.max_k_localTimer.startTimer();

            std::pair<std::vector<unsigned int>, ssize_t> newSeeds = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, localSolutionSpace
            );

            workingCandidateSet.first = newSeeds.first;
            workingCandidateSet.second = newSeeds.second;

            this->timeAggregator.max_k_localTimer.endTimer();
        }

        // spdlog::get("console")->info("all gather...");
        this->timeAggregator.allGatherTimer.startTimer();

        std::vector<unsigned int> globalAggregation;
        size_t totalData = this->cEngine.GatherPartialSolutions(workingCandidateSet, localSolutionSpace, globalAggregation, this->groupWorld);

        this->timeAggregator.allGatherTimer.endTimer();

        std::pair<std::vector<unsigned int>, size_t> globalCandidateSet;
        std::map<int, std::vector<int>> globalSolutionSpace;
        std::vector<std::pair<unsigned int, std::vector<unsigned int>>> candidateSets;

        // TODO: use group rank
        int groupRank;
        MPI_Comm_rank(this->groupWorld, &groupRank);
        if (groupRank == 0) {

            // spdlog::get("console")->info("unpacking local seeds in global process...");

            this->timeAggregator.allGatherTimer.startTimer();
            candidateSets = this->cEngine.aggregateLocalKSeeds(globalSolutionSpace, globalAggregation.data(), totalData);

            this->timeAggregator.allGatherTimer.endTimer();

            // spdlog::get("console")->info("global max_k_cover...");
            this->timeAggregator.max_k_globalTimer.startTimer();

            globalCandidateSet = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, globalSolutionSpace
            );

            candidateSets.push_back(std::make_pair(
                static_cast<unsigned int>(globalCandidateSet.second),
                globalCandidateSet.first
            ));

            this->timeAggregator.max_k_globalTimer.endTimer();
        }

        return {
            globalSolutionSpace,
            candidateSets
        };
    }
};