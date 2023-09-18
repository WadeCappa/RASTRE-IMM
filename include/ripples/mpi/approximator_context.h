#include "ripples/mpi/approximator_group.h"

template <typename GraphTy, typename ConfTy>
class ApproximatorContext {
    private:
    const std::vector<LazyLazyApproximatorGroup<GraphTy, ConfTy>*> groups;

    public:
    ApproximatorContext(const std::vector<LazyLazyApproximatorGroup<GraphTy, ConfTy>*> groups) : groups(groups) {}

    std::pair<std::vector<unsigned int>, unsigned int> getBestSeeds(
        const std::map<int, std::vector<int>> &localSolutionSpace,
        const unsigned int kprime,
        const size_t theta
    ) const {
        SolutionState state;
        state.solutionSpace = localSolutionSpace;

        for (size_t i = 0; i < this->groups.size(); i++) {
            std::cout << "working on level " << i << " which has solution space of size " << state.solutionSpace.size();
            std::cout << " and local solution utility of " << state.bestSolution.second;
	    std::cout << " and size of " << state.bestSolution.first.size() << std::endl;
            state = this->groups[i]->approximate(state, kprime, theta);
        }

        // std::cout << "returning best approximation" << std::endl;

        return state.bestSolution;
    }
};
