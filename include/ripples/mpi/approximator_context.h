#include "ripples/mpi/approximator_group.h"

class ApproximatorContext {
    private:
    const std::vector<ApproximatorGroups> &groups;

    void loadGroupResultAsAggregateSet(
        const std::pair<std::vector<unsigned int>, size_t>& groupResult, 
        std::map<int, std::vector<int>>& aggregateSets 
    ) {
        for (const auto & e : groupResult.first) {

        }
    }

    std::pair<std::vector<unsigned int>, unsigned int> getBestCandidate(
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

        return std::make_pair(bestCandidate, bestUtility);
    }

    public:

    ApproximatorContext(const std::vector<ApproximatorGroups> &groups) : groups(groups) {}

    std::pair<std::vector<unsigned int>, unsigned int> getBestSeeds(
        const std::map<int, std::vector<int>>& aggregateSets,
        const unsigned int kprime,
        const size_t theta
    ) const {
        // here is the flow; 
            // from each group object in order; 
                // allgather first
                // then you run seed selection. 
            // each group returns their best seeds. 
            // if you are not the group's leader, then you are on the last group in this->groups, you are done looping, exit
            // otherwise you have data returned by group->approximate, feed this into the next group
        
        std::map<int, std::vector<int>> currentData = aggregateSets; 
        std::pair<std::vector<unsigned int>, unsigned int> localCandidateSet;
        for (const auto & group : groups) {
            SolutionCandidateSets candidates = group.approximate(currentData, localCandidateSet, kprime, theta);
            localCandidateSet = getBestCandidate(candidates.localCandidateSets); 
            currentData = candidates.allKMSeeds;
        }

        return localCandidateSet; 
    }
};