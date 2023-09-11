

#include "ripples/generate_rrr_sets.h"

class MartingleBuilder {
    public:
    template <
        typename GraphTy,
        typename ConfTy>
    static std::vector<MPI_Comm> buildCommGroups(
        const unsigned int levels,
        const unsigned int branchingFactor,
        const int rank,
        const int worldSize
    ) {
        std::vector<MPI_Comm> groups;
        int currentRank = rank;
        unsigned int currentLevel = 0;
        unsigned int currentWorldSize = worldSize;
        bool addingNewGroups = true;

        for (unsigned int level = 0; level < levels; level++) {
            unsigned int numberOfGroups = std::ceil(currentWorldSize / branchingFactor);

            MPI_Comm newComm;
            // TODO: This creates lobsided groups when branchingFactor is not a factor of m^(1/levels)
            int color = std::floor(currentRank / numberOfGroups);
            MPI_Comm_split(
                MPI_COMM_WORLD, 
                addingNewGroups ? color : MPI_UNDEFINED,
                rank, 
                &newComm
            );

            currentWorldSize = currentWorldSize / branchingFactor;

            if (addingNewGroups) {
                groups.push_back(newComm);
                int groupRank;
                MPI_Comm_rank(newComm, &groupRank);
                if (groupRank == 0) {
                    currentRank = currentRank / branchingFactor;
                } else {
                    addingNewGroups = false;
                    currentRank = -1;
                }
            }
        }

        return groups;
    }

    template <
        typename GraphTy,
        typename ConfTy>
    static void buildApproximatorGroups(
        std::vector<LazyLazyApproximatorGroup<GraphTy, ConfTy>*> &approximatorGroups, // modified
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

    // template <
    //     typename GraphTy,
    //     typename ConfTy,
    //     typename diff_model_tag,
    //     typename RRRGeneratorTy>
    // static MartingaleContext<GraphTy, ConfTy, RRRGeneratorTy> build(
    //     const GraphTy &G, 
    //     const ConfTy &CFG, 
    //     double l_value, 
    //     RRRGeneratorTy &gen,
    //     ripples::IMMExecutionRecord &record, 
    //     diff_model_tag &model_tag, 
    //     unsigned int levels,
    //     unsigned int branchingFactor) 
    // {
    //     // no sequential version available
    //     using execution_tag = ripples::omp_parallel_tag;
    //     using vertex_type = typename GraphTy::vertex_type;

    //     CommunicationEngine<GraphTy> cEngine = CommunicationEngineBuilder<GraphTy>::BuildCommunicationEngine();
    //     TransposeRRRSets<GraphTy> tRRRSets(G.num_nodes());
    //     TimerAggregator timeAggregator;
    //     std::vector<int> vertexToProcess(cEngine.DistributeVertices(CFG.use_streaming, G));

    //     DefaultSampler<GraphTy, diff_model_tag, RRRGeneratorTy, execution_tag> sampler(
    //         cEngine.GetSize(), G, gen, record, model_tag
    //     );

    //     OwnershipManager<GraphTy> ownershipManager(G.num_nodes(), cEngine, vertexToProcess);

    //     std::vector<MPI_Comm> groups = buildCommGroups<GraphTy, ConfTy>(
    //         levels, branchingFactor, cEngine.GetRank(), cEngine.GetSize()
    //     );

    //     std::vector<LazyLazyApproximatorGroup<GraphTy, ConfTy>*> approximatorGroups;
    //     buildApproximatorGroups(approximatorGroups, groups, CFG, vertexToProcess, timeAggregator, cEngine);

    //     ApproximatorContext<GraphTy, ConfTy> approximator(approximatorGroups);

    //     MartingaleContext<GraphTy, ConfTy, RRRGeneratorTy, diff_model_tag, execution_tag> martingaleContext(
    //         sampler, ownershipManager, approximator, G, CFG, l_value, record, cEngine, timeAggregator
    //     );

    //     return martingaleContext;
    // }
};