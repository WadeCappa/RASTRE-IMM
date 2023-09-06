#include "ripples/mpi/approximator_group.h"

class ApproximatorContext {
    private:
    const std::vector<ApproximatorGroup> &groups;

    std::pair<std::vector<unsigned int>, ssize_t> getBestCandidate(
        const std::vector<std::pair<unsigned int, std::vector<unsigned int>>> &allCandidateSets
    ) const {
        ssize_t bestUtility = 0;
        std::vector<unsigned int> bestCandidate;

        for (const auto & e : allCandidateSets) {
            if (e.first > bestUtility) {
                bestCandidate = e.second;
                bestUtility = e.first;
            }
        }

        return std::make_pair(bestCandidate, bestUtility);
    }

    public:
    ApproximatorContext(const std::vector<ApproximatorGroup> &groups) : groups(groups) {}

    std::pair<std::vector<unsigned int>, unsigned int> getBestSeeds(
        const std::map<int, std::vector<int>>& localSolutionSpace,
        const unsigned int kprime,
        const size_t theta
    ) const {
        std::map<int, std::vector<int>> currentData = localSolutionSpace; 
        std::pair<std::vector<unsigned int>, ssize_t> localCandidateSet;
        for (const auto & group : groups) {
            SolutionCandidateSets candidates = group.approximate(currentData, localCandidateSet, kprime, theta);
            localCandidateSet = getBestCandidate(candidates.localCandidateSets); 
            currentData = candidates.allKMSeeds;
        }

        return localCandidateSet; 
    }
};