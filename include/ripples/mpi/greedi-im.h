//===------------------------------------------------------------*- C++ -*-===//
//
//             Ripples: A C++ Library for Influence Maximization
//                  Marco Minutoli <marco.minutoli@pnnl.gov>
//                   Pacific Northwest National Laboratory
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2019, Battelle Memorial Institute
//
// Battelle Memorial Institute (hereinafter Battelle) hereby grants permission
// to any person or entity lawfully obtaining a copy of this software and
// associated documentation files (hereinafter “the Software”) to redistribute
// and use the Software in source and binary forms, with or without
// modification.  Such person or entity may use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and may permit
// others to do so, subject to the following conditions:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimers.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Other than as used herein, neither the name Battelle Memorial Institute or
//    Battelle may be used in any form whatsoever without the express written
//    consent of Battelle.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//===----------------------------------------------------------------------===//

#ifndef RIPPLES_MPI_IMM_H
#define RIPPLES_MPI_IMM_H

#include "mpi.h"

#include <cstddef>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <thread>

#include "trng/lcg64.hpp"

#include "ripples/generate_rrr_sets.h"
#include "ripples/imm.h"
#include "ripples/imm_execution_record.h"
#include "ripples/mpi/find_most_influential.h"
#include "ripples/utility.h"

#include "ripples/generate_rrr_sets.h"
#include "ripples/mpi/communication_engine.h"
#include "ripples/TimerAggregator.h"
#include "ripples/max_k_cover.h"

#include <time.h>
#include <cstdlib>
#include <numeric>
#include <random>
#include "ripples/imm_execution_record.h"

#include "ripples/mpi/streaming_engine.h"
#include "ripples/mpi/approximator_context.h"
#include "ripples/mpi/ownership_manager.h"
#include "ripples/sampler_context.h"
#include "ripples/mpi/martingale_context.h"
#include "ripples/mpi/martingale_builder.h"

namespace ripples {
namespace mpi {

template <typename ex_tag>
struct MPI_Plus_X {
  // using generate_ex_tag
  // using seed_selection_ex_tag
};

template <>
struct MPI_Plus_X<mpi_omp_parallel_tag> {
  using generate_ex_tag = omp_parallel_tag;
  using seed_selection_ex_tag = mpi_omp_parallel_tag;
};

//! Split a random number generator into one sequence per MPI rank.
//!
//! \tparam PRNG The type of the random number generator.
//!
//! \param gen The parallel random number generator to split.
template <typename PRNG>
void split_generator(PRNG &gen) {
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  gen.split(world_size, world_rank);
}

template <typename PRNG>
std::vector<PRNG> rank_split_generator(const PRNG &gen) {
  size_t max_num_threads(1);

#pragma omp single
  max_num_threads = omp_get_max_threads();

  std::vector<trng::lcg64> generator(max_num_threads, gen);

#pragma omp parallel
  {
    generator[omp_get_thread_num()].split(omp_get_num_threads(),
                                          omp_get_thread_num());
  }
  return generator;
}


template <typename GraphTy, typename ConfTy, typename diff_model_tag,
          typename RRRGeneratorTy, typename ExTagTrait>
auto run_greedimm(
  const GraphTy &G, 
  const ConfTy &CFG, 
  TimerAggregator &timeAggregator,
  double l_value, 
  RRRGeneratorTy &gen,
  IMMExecutionRecord &record, 
  diff_model_tag &&model_tag, 
  ExTagTrait &&
) 
{
  // no sequential version available
  using execution_tag = ripples::omp_parallel_tag;
  using vertex_type = typename GraphTy::vertex_type;

  CommunicationEngine<GraphTy> cEngine = CommunicationEngineBuilder<GraphTy>::BuildCommunicationEngine(CFG.use_streaming);
  TransposeRRRSets<GraphTy> tRRRSets(G.num_nodes());
  std::vector<int> vertexToProcess(cEngine.DistributeVertices(CFG.use_streaming, G));

  DefaultSampler<GraphTy, diff_model_tag, RRRGeneratorTy, execution_tag> sampler(
      cEngine.GetSize(), G, gen, record, model_tag
  );

  OwnershipManager<GraphTy> ownershipManager(G.num_nodes(), cEngine, vertexToProcess);

  std::vector<ApproximatorGroup*> approximatorGroups;
  approximatorGroups.push_back(new StreamingApproximatorGroup<GraphTy, ConfTy>(MPI_COMM_WORLD, vertexToProcess, CFG, cEngine));
  std::vector<ApproximatorContext> approximators;
  approximators.push_back(ApproximatorContext (approximatorGroups, timeAggregator));

  MartingaleContext<GraphTy, ConfTy, RRRGeneratorTy, diff_model_tag, execution_tag> martingaleContext(
      sampler, ownershipManager, approximators, G, CFG, l_value, record, cEngine, timeAggregator
  );

  auto res = martingaleContext.approximateInfMax();

  return res;
}

template <typename GraphTy, typename ConfTy, typename diff_model_tag,
          typename RRRGeneratorTy, typename ExTagTrait>
auto run_randgreedi(
  const GraphTy &G, 
  const ConfTy &CFG, 
  TimerAggregator &timeAggregator,
  double l_value, 
  RRRGeneratorTy &gen,
  IMMExecutionRecord &record, 
  diff_model_tag &&model_tag, 
  const std::vector<unsigned int> &branchingFactors,
  ExTagTrait &&
) 
{
  // no sequential version available
  using execution_tag = ripples::omp_parallel_tag;
  using vertex_type = typename GraphTy::vertex_type;

  CommunicationEngine<GraphTy> cEngine = CommunicationEngineBuilder<GraphTy>::BuildCommunicationEngine(CFG.use_streaming);
  TransposeRRRSets<GraphTy> tRRRSets(G.num_nodes());
  std::vector<int> vertexToProcess(cEngine.DistributeVertices(CFG.use_streaming, G));

  DefaultSampler<GraphTy, diff_model_tag, RRRGeneratorTy, execution_tag> sampler(
      cEngine.GetSize(), G, gen, record, model_tag
  );

  OwnershipManager<GraphTy> ownershipManager(G.num_nodes(), cEngine, vertexToProcess);

  std::vector<ApproximatorContext> approximators = MartingleBuilder::buildApproximatorContexts<GraphTy, ConfTy>(branchingFactors, cEngine.GetSize(), CFG, timeAggregator, vertexToProcess, cEngine);

  MartingaleContext<GraphTy, ConfTy, RRRGeneratorTy, diff_model_tag, execution_tag> martingaleContext(
      sampler, ownershipManager, approximators, G, CFG, l_value, record, cEngine, timeAggregator
  );

  auto res = martingaleContext.approximateInfMax();

  return res;
}



}  // namespace mpi
}  // namespace ripples

#endif  // RIPPLES_MPI_IMM_H
