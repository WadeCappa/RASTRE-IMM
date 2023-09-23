typedef struct solutionState {
    const std::map<int, std::vector<int>> &solutionSpace;
    const std::pair<std::vector<unsigned int>, unsigned int> &bestSolution;
} SolutionState ;

#include "ripples/mpi/approximator_group.h"


class ApproximatorContext {
    private:
    const std::vector<ApproximatorGroup*> groups;

    public:
    ApproximatorContext(const std::vector<ApproximatorGroup*> groups) : groups(groups) {}

    std::pair<std::vector<unsigned int>, unsigned int> getBestSeeds(
        const std::map<int, std::vector<int>> &initialSolutionSpace,
        const unsigned int kprime,
        const size_t theta
    ) const {
        SolutionState initialState = {initialSolutionSpace, std::make_pair(std::vector<unsigned int>(), 0)};
        const SolutionState *currentState = &initialState;
        std::vector<std::map<int, std::vector<int>>> solutionSpaces(this->groups.size());
        std::vector<std::pair<std::vector<unsigned int>, unsigned int>> bestSolutions(this->groups.size());
        std::vector<SolutionState> solutionStates;

        for (size_t i = 0; i < this->groups.size(); i++) {
            auto nextSolutionSpace = &(solutionSpaces[i]);
            auto nextBestSolution = &(bestSolutions[i]);

            // std::cout << "working on level " << i << " which has solution space of size " << currentState->solutionSpace.size();
            // std::cout << " and local solution utility of " << currentState->bestSolution.second;
            // std::cout << " and size of " << currentState->bestSolution.first.size() << std::endl;

            this->groups[i]->approximate(*currentState, *nextSolutionSpace, *nextBestSolution, kprime, theta);
            solutionStates.push_back({*nextSolutionSpace, *nextBestSolution});
            currentState = &(solutionStates[i]);
        }

        return currentState->bestSolution;
    }
};
