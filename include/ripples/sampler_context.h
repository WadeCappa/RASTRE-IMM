#include "ripples/generate_rrr_sets.h"
#include "ripples/imm_execution_record.h"
#include "ripples/mpi/find_most_influential.h"
#include "ripples/utility.h"
#include "ripples/imm_execution_record.h"
#include "ripples/mpi/imm.h"

// template <typename GraphTy>
// class SamplerContext {
//     private:

//     public:
//     virtual void addNewSamples(
//         TransposeRRRSets<GraphTy> &tRRRSets, 
//         const size_t counterStart, 
//         const size_t delta
//     ) = 0;
// };

template <
    typename GraphTy,
    typename diff_model_tag,
    typename RRRGeneratorTy,
    typename execution_tag
>
class DefaultSampler {
    private:
    ripples::omp_parallel_tag ex_tag;
    diff_model_tag &model_tag;

    const unsigned int world_size;
    const GraphTy &graph;
    RRRGeneratorTy &gen;
    ripples::IMMExecutionRecord &record;

    public:
    DefaultSampler(
        unsigned int world_size,
        const GraphTy &input_G, 
        RRRGeneratorTy &gen,
        ripples::IMMExecutionRecord &record,
        diff_model_tag &model_tag
    ) : world_size(world_size), graph(input_G), gen(gen), record(record), model_tag(model_tag) {
    }

    // Non-const, but because we want random generation this is fine.
    void addNewSamples(
        TransposeRRRSets<GraphTy> &tRRRSets, 
        const size_t counterStart, 
        const size_t delta
    ) {
        GenerateTransposeRRRSets(
            tRRRSets, counterStart, delta, 
            this->graph, this->gen, this->record,
            std::forward<diff_model_tag>(this->model_tag),
            std::forward<execution_tag>(this->ex_tag)
        );
    }
};