#include "ripples/mpi/approximator_context.h"
#include "ripples/mpi/seed_set_aggregator.h"
#include "ripples/sampler_context.h"

class MargingaleRunner {
    private:
    SamplerContext &sampler;
    const ApproximatorContext &approximator;
    const SeedSetAggregator &aggregator;

    const GraphTy &G;
    const ConfTy &CFG;
    diff_model_tag &model_tag;
    RRRGeneratorTy &gen;
    IMMExecutionRecord &record;
    const CommunicationEngine<GraphTy> &cEngine;

    const double l;
    size_t previous_theta = 0;
    size_t thetaPrime = 0;

    std::map<int, std::vector<int>> aggregateSets;
    std::vector<size_t> old_sampling_sizes;
    RRRsetAllocator<typename GraphTy::vertex_type> allocator;
    TransposeRRRSets<GraphTy> tRRRSets;
    const std::vector<int> vertexToProcess;

    TimerAggregator &timeAggregator;

    ssize_t getThetaPrime(
        ssize_t x, double epsilonPrime, double l, size_t k, size_t num_nodes 
    ) {
        k = std::min(k, num_nodes/2);
        return (2 + 2. / 3. * epsilonPrime) *
            (l * std::log(num_nodes) + logBinomial(num_nodes, k) +
            std::log(std::log2(num_nodes))) *
            std::pow(2.0, x) / (epsilonPrime * epsilonPrime);
    }

    public:
    MargingaleRunner(
        SamplerContext &sampler,
        const SeedSetAggregator &aggregator,
        const ApproximatorContext &approximator,

        const GraphTy &input_G, 
        const ConfTy &input_CFG,
        ripples::omp_parallel_tag &e_tag,
        diff_model_tag &model_tag,
        const double input_l,
        RRRGeneratorTy &gen,
        IMMExecutionRecord &record,
        const CommunicationEngine<GraphTy> &cEngine,
        TimerAggregator &timeAggregator
    ) 
        : 
            approximator(approximator), 
            aggregator(aggregator), 
            sampler(sampler),

            G(input_G),
            CFG(input_CFG), 
            ex_tag(e_tag), 
            model_tag(model_tag),
            tRRRSets(G.num_nodes()), 
            l(input_l * (1 + 1 / std::log2(G.num_nodes()))),
            vertexToProcess(cEngine.DistributeVertices(input_CFG.use_streaming, input_G)),
            old_sampling_sizes(input_G.num_nodes(), 0),
            gen(gen),
            record(record),
            cEngine(cEngine),
            timeAggregator(timeAggregator)
        {
    }
}