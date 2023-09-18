
typedef struct solutionState {
    std::map<int, std::vector<int>> solutionSpace;
    std::pair<std::vector<unsigned int>, unsigned int> bestSolution;
} SolutionState ;

class ApproximatorGroup {
    protected: 
    const MPI_Comm groupWorld;
    const std::vector<int> &vertexToProcess;

    template <typename GraphTy>
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
		
	// std::cout << "getting best solution from " << allCandidateSets.size() << " options" << std::endl;

        for (const auto & e : allCandidateSets) {
            if (e.first > bestUtility) {
                bestCandidate = e.second;
                bestUtility = e.first;
            }
        }

        return std::make_pair(bestCandidate, bestUtility);
    }

    public:
    ApproximatorGroup(
        MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess
    ) : vertexToProcess(vertexToProcess), groupWorld(groupWorld) {}

    virtual SolutionState approximate(
        const SolutionState &previousState,
        const int kprime, 
        const size_t theta
    ) = 0;
};

template <typename GraphTy, typename ConfTy>
class StreamingApproximatorGroup : public ApproximatorGroup { 
    private:
    TimerAggregator &timeAggregator;
    const ConfTy &CFG;
    const CommunicationEngine<GraphTy> &cEngine;

    public:
    StreamingApproximatorGroup(
        MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess,
        TimerAggregator &timeAggregator,
        const ConfTy &input_CFG,
        const CommunicationEngine<GraphTy> &cEngine
    ) : ApproximatorGroup(groupWorld, vertexToProcess), timeAggregator(timeAggregator), CFG(input_CFG), cEngine(cEngine) {}

    SolutionState approximate(
        const SolutionState &previousState,
        const int kprime, 
        const size_t theta
    ) override {
        SolutionState currentState;
        int worldSize;
        int groupRank;
        MPI_Comm_rank(this->groupWorld, &groupRank);
        MPI_Comm_size(this->groupWorld, &worldSize);
        if (groupRank == 0) {
            std::pair<std::vector<unsigned int>, size_t> globalCandidateSet;
            std::map<int, std::vector<int>> newSolutionSpace;
            StreamingRandGreedIEngine<GraphTy> streamingEngine(this->CFG.k, kprime, theta, (double)this->CFG.epsilon_2, worldSize - 1, this->cEngine, this->timeAggregator);
            globalCandidateSet = streamingEngine.Stream(newSolutionSpace);

            currentState.bestSolution = globalCandidateSet;
            currentState.solutionSpace = newSolutionSpace;
        } else {
            this->timeAggregator.totalSendTimer.startTimer();
            
            StreamingMaxKCover<GraphTy> localKCoverEngine(this->CFG.k, kprime, theta, this->timeAggregator, this->cEngine);
            localKCoverEngine.useLazyGreedy();
            localKCoverEngine.run_max_k_cover(previousState.solutionSpace);

            this->timeAggregator.totalSendTimer.endTimer();
        }

        return currentState;
    }
};

template <typename GraphTy, typename ConfTy>
class LazyLazyApproximatorGroup : public ApproximatorGroup {
    private:
    TimerAggregator &timeAggregator;
    const ConfTy &CFG;
    const CommunicationEngine<GraphTy> &cEngine;

    public:
    LazyLazyApproximatorGroup(
        MPI_Comm groupWorld,
        const std::vector<int> &vertexToProcess,
        TimerAggregator &timeAggregator,
        const ConfTy &input_CFG,
        const CommunicationEngine<GraphTy> &cEngine
    ) : ApproximatorGroup(groupWorld, vertexToProcess), timeAggregator(timeAggregator), CFG(input_CFG), cEngine(cEngine) {}

    SolutionState approximate(
        const SolutionState &previousState,
        const int kprime, 
        const size_t theta
    ) override {
        SolutionState currentState;
        currentState.bestSolution.first = previousState.bestSolution.first;
        currentState.bestSolution.second = previousState.bestSolution.second;

        if (currentState.bestSolution.first.size() == 0) {
            // std::cout << "init of local solution using world size of " << previousState.solutionSpace.size() << std::endl;;
            this->timeAggregator.max_k_localTimer.startTimer();

            std::pair<std::vector<unsigned int>, unsigned int> newSeeds = this->SolveKCover<GraphTy>(
                this->CFG.k, kprime, theta, this->timeAggregator, previousState.solutionSpace
            );

            currentState.bestSolution.second = newSeeds.second;
            currentState.bestSolution.first = std::vector<unsigned int>(newSeeds.first.begin(), newSeeds.first.end());

            this->timeAggregator.max_k_localTimer.endTimer();
        }

        // std::cout << "about to do gather on " << currentState.bestSolution.first.size();
        // std::cout << " seeds of utility " << currentState.bestSolution.second;
        // std::cout << " using solution space of size " << previousState.solutionSpace.size();
        // std::cout << " for rank " << cEngine.GetRank() << std::endl;

        // spdlog::get("console")->info("all gather...");
        this->timeAggregator.allGatherTimer.startTimer();

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

        int globalRank; 
        MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

        int groupRank;
        MPI_Comm_rank(this->groupWorld, &groupRank);
        if (groupRank == 0) {

            std::cout << "global rank " << globalRank << " is group leader" << std::endl;

            this->timeAggregator.allGatherTimer.startTimer();
            candidateSets = this->cEngine.aggregateLocalKSeeds(globalSolutionSpace, globalAggregation.data(), totalData);

            this->timeAggregator.allGatherTimer.endTimer();

            // spdlog::get("console")->info("global max_k_cover...");
            this->timeAggregator.max_k_globalTimer.startTimer();

            // std::cout << "global process getting solution from " << globalSolutionSpace.size() << " possible solutions " << std::endl;

            globalCandidateSet = this->SolveKCover<GraphTy>(
                this->CFG.k, kprime, theta, this->timeAggregator, globalSolutionSpace
            );

            candidateSets.push_back(std::make_pair(
                globalCandidateSet.second,
                globalCandidateSet.first
            ));

            // currentState.bestSolution.first = globalCandidateSet.first;
            // currentState.bestSolution.second = globalCandidateSet.second;

            this->timeAggregator.max_k_globalTimer.endTimer();
        }

        currentState.solutionSpace = globalSolutionSpace;
        currentState.bestSolution = this->getBestCandidate(candidateSets);
        // std::cout << "got best candidate of " << currentState.bestSolution.second << " and size " << currentState.bestSolution.first.size() << " for rank " << globalRank << std::endl;

        return currentState;
    }
};
