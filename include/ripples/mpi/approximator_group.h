
typedef struct solutionState {
    std::map<int, std::vector<int>> solutionSpace;
    std::pair<std::vector<unsigned int>, unsigned int> bestSolution;
} SolutionState ;

class ApproximatorGroup {
    protected: 
    const MPI_Comm groupWorld;
    const std::vector<int> &vertexToProcess;

    public:
    ApproximatorGroup(
        const MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess
    ) : vertexToProcess(vertexToProcess), groupWorld(groupWorld) {}

    virtual void approximate(
        SolutionState &solutionSets, // modified as side effect
        const std::pair<std::vector<unsigned int>, unsigned int> &previousGroupSolution,
        const int kprime, 
        const size_t theta
    ) = 0;
};

// class DummyApproximatorGroup : public ApproximatorGroup {
//     public:
//     DummyApproximatorGroup(
//         const MPI_Comm groupWorld,
//         const std::vector<int> &vertexToProcess
//     ) : ApproximatorGroup(groupWorld, vertexToProcess) {}

//     void approximate(
//         SolutionState &solutionSets, // modified as side effect
//         const std::pair<std::vector<unsigned int>, unsigned int> &previousGroupSolution,
//         const int kprime, 
//         const size_t theta
//     ) {
//         std::cout << "inside dummy approximator" << std::endl;

//         int group_rank; 
//         MPI_Comm_rank(groupWorld, &group_rank);
        
//         if (group_rank == 0) {
//             solutionSets.solutionSpace = std::map<int, std::vector<int>>{
//                     {0, std::vector<int>{0,1,2,3,4,5,6,7,8,9}},
//                     {1, std::vector<int>{0,1,2,3,4,5,6,7,8,9}},
//                     {2, std::vector<int>{0,1,2,3,4,5,6,7,8,9}}
//                 };
//             solutionSets.localCandidateSets = std::vector<std::pair<unsigned int, std::vector<unsigned int>>>{
//                     std::make_pair(80000, std::vector<unsigned int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9})
//                 };
//         } else {
//             solutionSets.solutionSpace = std::map<int, std::vector<int>>();
//             solutionSets.localCandidateSets = std::vector<std::pair<unsigned int, std::vector<unsigned int>>>();
//         }
//     }
// };

template <typename GraphTy, typename ConfTy>
class LazyLazyApproximatorGroup {
    private:
    TimerAggregator &timeAggregator;
    const ConfTy &CFG;
    const MPI_Comm groupWorld;
    const std::vector<int> &vertexToProcess;
    const CommunicationEngine<GraphTy> &cEngine;

    static std::pair<std::vector<unsigned int>, unsigned int> SolveKCover(
        const int k,
        const int kprime,
        const size_t max_element, 
        TimerAggregator &timeAggregator,
        const std::map<int, std::vector<int>>& elements
    ) {
        MaxKCover<GraphTy> localKCoverEngine(k, kprime, max_element, timeAggregator);
        localKCoverEngine.useLazyGreedy();
        std::pair<std::vector<unsigned int>, size_t> temp = localKCoverEngine.run_max_k_cover(elements);
        // std::cout << "got local seeds" << std::endl;
        return std::make_pair(
            temp.first,
            static_cast<int>(temp.second)
        );
    }

    static std::pair<std::vector<unsigned int>, unsigned int> getBestCandidate(
        const std::vector<std::pair<unsigned int, std::vector<unsigned int>>> &allCandidateSets
    ) {
        unsigned int bestUtility = 0;
        std::vector<unsigned int> bestCandidate;

        for (const auto & e : allCandidateSets) {
            if (e.first > bestUtility) {
                bestCandidate = e.second;
                bestUtility = e.first;
            }
        }

        std::cout << "got best candidate of " << bestUtility << " and size " << bestCandidate.size() << std::endl;

        return std::make_pair(bestCandidate, bestUtility);
    }

    public:
    LazyLazyApproximatorGroup(
        MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess,
        TimerAggregator &timeAggregator,
        const ConfTy &input_CFG,
        const CommunicationEngine<GraphTy> &cEngine
    ) : timeAggregator(timeAggregator), CFG(input_CFG), cEngine(cEngine), groupWorld(groupWorld), vertexToProcess(vertexToProcess) {}

    SolutionState approximate(
        const SolutionState &previousState,
        const int kprime, 
        const size_t theta
    ) {
        SolutionState currentState;
        currentState.bestSolution = previousState.bestSolution;

        if (currentState.bestSolution.first.size() == 0) {
            std::cout << "init of local solution" << std::endl;
            this->timeAggregator.max_k_localTimer.startTimer();

            std::pair<std::vector<unsigned int>, unsigned int> newSeeds = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, previousState.solutionSpace
            );

            // TODO: Fix this, should not be assigning first to second and second to first
            currentState.bestSolution.second = newSeeds.second;
            currentState.bestSolution.first = newSeeds.first;

            this->timeAggregator.max_k_localTimer.endTimer();
        }

        // spdlog::get("console")->info("all gather...");
        this->timeAggregator.allGatherTimer.startTimer();

        std::cout << "aggregating current best solution of " << currentState.bestSolution.first.size() << " and utility of " << currentState.bestSolution.second << std::endl;
        std::vector<unsigned int> globalAggregation;
        size_t totalData = this->cEngine.GatherPartialSolutions(
            currentState.bestSolution,
            previousState.solutionSpace, 
            globalAggregation, 
            this->groupWorld
        );

        this->timeAggregator.allGatherTimer.endTimer();

        std::pair<std::vector<unsigned int>, size_t> globalCandidateSet;
        std::map<int, std::vector<int>> globalSolutionSpace;
        std::vector<std::pair<unsigned int, std::vector<unsigned int>>> candidateSets;

        // TODO: use group rank
        int groupRank;
        MPI_Comm_rank(this->groupWorld, &groupRank);
        if (groupRank == 0) {

            int globalRank; 
            MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

            std::cout << "global rank " << globalRank << " is group leader" << std::endl;

            this->timeAggregator.allGatherTimer.startTimer();
            candidateSets = this->cEngine.aggregateLocalKSeeds(globalSolutionSpace, globalAggregation.data(), totalData);

            this->timeAggregator.allGatherTimer.endTimer();

            // spdlog::get("console")->info("global max_k_cover...");
            this->timeAggregator.max_k_globalTimer.startTimer();

            std::cout << "solving for group globally" << std::endl;
            globalCandidateSet = this->SolveKCover(
                this->CFG.k, kprime, theta, this->timeAggregator, globalSolutionSpace
            );

            candidateSets.push_back(std::make_pair(
                static_cast<unsigned int>(globalCandidateSet.second),
                globalCandidateSet.first
            ));

            this->timeAggregator.max_k_globalTimer.endTimer();
        }

        std::cout << "rank " << this->cEngine.GetRank() << " has returned " << candidateSets.size() << " candidate sets, and a solution space of " << globalSolutionSpace.size() << std::endl;

        currentState.solutionSpace = globalSolutionSpace;
        currentState.bestSolution = this->getBestCandidate(candidateSets);
        return currentState;
    }
};