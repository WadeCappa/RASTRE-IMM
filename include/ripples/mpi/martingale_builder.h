

#include "ripples/generate_rrr_sets.h"

class MartingleBuilder {
    public:

    template <
        typename GraphTy,
        typename ConfTy>
    static ApproximatorContext buildApproximatorContext(
        const unsigned int levels,
        const unsigned int branchingFactor,
        const ConfTy &CFG, 
        TimerAggregator &timeAggregator,
        const std::vector<int> &vertexToProcess,
        const CommunicationEngine<GraphTy> &cEngine
    ) {
        std::vector<MPI_Comm> groups = MartingleBuilder::buildCommGroups<GraphTy, ConfTy>(
            levels, branchingFactor, cEngine.GetRank(), cEngine.GetSize()
        );

        std::vector<ApproximatorGroup*> approximatorGroups;
        MartingleBuilder::buildApproximatorGroups(approximatorGroups, groups, CFG, vertexToProcess, timeAggregator, cEngine);

        return ApproximatorContext(approximatorGroups);
    }

    template <
        typename GraphTy,
        typename ConfTy>
    static std::vector<MPI_Comm> buildCommGroups(
        const unsigned int levels, 
        const unsigned int branchingFactor, // 4
        const int rank,
        const int worldSize // 8
    ) {
        std::vector<MPI_Comm> groups;
        unsigned int currentWorldSize = worldSize;
        bool addingNewGroups = true; 

        for (unsigned int level = 0; level < levels; level++) {
            // std::cout << "starting level " << level << " at rank " << rank << std::endl;
            unsigned int numberOfGroups = std::ceil((double)currentWorldSize / (double)branchingFactor); 

            if (numberOfGroups == 0) {
                std::cout << "ROUNDING ERROR, SHOULD NEVER SEE THIS STATMENT" << std::endl;
                exit(1);
            }

            MPI_Comm newComm;
            // TODO: This creates lobsided groups when branchingFactor is not a factor of m^(1/levels)
            int color = rank % numberOfGroups;
            // std::cout << "level " << level << " has number of groups " << numberOfGroups;
            // std::cout << " and rank " << rank << " has color " << color << std::endl;
            MPI_Comm_split(
                MPI_COMM_WORLD, 
                addingNewGroups ? color : MPI_UNDEFINED,
                rank, 
                &newComm
            );

            currentWorldSize = numberOfGroups;

            if (addingNewGroups) {
                groups.push_back(newComm);

                int groupRank;
                MPI_Comm_rank(newComm, &groupRank);
                if (groupRank != 0) {
                    addingNewGroups = false;
                }
            }
        }

        return groups;
    }

    template <
        typename GraphTy,
        typename ConfTy>
    static void buildApproximatorGroups(
        std::vector<ApproximatorGroup*> &approximatorGroups, // modified
        const std::vector<MPI_Comm> &commGroups,
        const ConfTy &CFG, 
        const std::vector<int> &vertexToProcess,
        TimerAggregator &timeAggregator,
        const CommunicationEngine<GraphTy> &cEngine) {
        for (const auto e : commGroups) {
            // ApproximatorGroup* newGroup = new DummyApproximatorGroup(e, vertexToProcess);

            approximatorGroups.push_back(new LazyLazyApproximatorGroup<GraphTy, ConfTy>(
                e, vertexToProcess, timeAggregator, CFG, cEngine
            ));
        }
    }
};
