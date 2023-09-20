typedef struct solutionState {
    std::map<int, std::vector<int>> solutionSpace;
    std::pair<std::vector<unsigned int>, unsigned int> bestSolution;
} SolutionState ;

#include "ripples/mpi/approximator_group.h"


class ApproximatorContext {
    private:
    const std::vector<ApproximatorGroup*> groups;

    public:
    ApproximatorContext(const std::vector<ApproximatorGroup*> groups) : groups(groups) {}

    void getBestSeeds(
        SolutionState &inputState,
        const unsigned int kprime,
        const size_t theta
    ) const {
        inputState.bestSolution = std::make_pair(std::vector<unsigned int>(), 0);
        const SolutionState *currentState = &inputState;
        std::vector<SolutionState> states(this->groups.size());

        for (size_t i = 0; i < this->groups.size(); i++) {
            SolutionState *nextState = &(states[i]);
            // std::cout << "working on level " << i << " which has solution space of size " << state.solutionSpace.size();
            // std::cout << " and local solution utility of " << state.bestSolution.second;
            // std::cout << " and size of " << state.bestSolution.first.size() << std::endl;

            this->groups[i]->approximate(*currentState, *nextState, kprime, theta);
            currentState = nextState;
        }

        inputState.bestSolution = currentState->bestSolution;
    }
};
