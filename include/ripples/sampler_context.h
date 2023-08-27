#include "ripples/generate_rrr_sets.h"

class SamplerContext {
    private:

    public:
    void addNewSamples(const size_t theta) = 0;
}

class DefaultSampler : public SamplerContext {
    private:
    ripples::omp_parallel_tag &ex_tag;
    const ripples::omp_parallel_tag e_tag = ripples::mpi::MPI_Plus_X<ripples::mpi_omp_parallel_tag>{};

    size_t previous_theta;
    const unsigned int world_size;
    TransposeRRRSets<GraphTy> &tRRRSets;
    const GraphTy &graph;
    RRRGeneratorTy &gen;
    IMMExecutionRecord &record;

    public:
    DefaultSampler(
        TransposeRRRSets<GraphTy> &tRRRSets,
        unsigned int world_size,
        const GraphTy &input_G, 
        RRRGeneratorTy &gen,
        IMMExecutionRecord &record,
        diff_model_tag &model_tag
    ) : tRRRSets(tRRRSets), world_size(world_size), graph(input_G), gen(gen), record(record) {
        this->previous_theta = 0;
    }

    void addNewSamples(const size_t theta) override {
        size_t delta = (theta - previous_theta) / this->world_size;

        GenerateTransposeRRRSets(
            this->tRRRSets, this->previous_theta, delta, 
            this->graph, this->gen, this->record,
            std::forward<diff_model_tag>(this->model_tag),
            std::forward<execution_tag>(this->ex_tag)
        );

        this->previous_theta = theta;
    }
}