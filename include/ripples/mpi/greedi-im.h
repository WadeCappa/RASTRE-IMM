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

#include "trng/lcg64.hpp"

#include "ripples/generate_rrr_sets.h"
#include "ripples/imm.h"
#include "ripples/imm_execution_record.h"
#include "ripples/mpi/find_most_influential.h"
#include "ripples/utility.h"

#include "ripples/generate_rrr_sets.h"
#include "ripples/mpi/communication_engine.h"
#include "ripples/max_k_cover.h"
#include "ripples/TimerAggregator.h"

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

//! Compute ThetaPrime for the MPI implementation.
//!
//! \param x The index of the current iteration.
//! \param epsilonPrime Parameter controlling the approximation factor.
//! \param l Parameter usually set to 1.
//! \param k The size of the seed set.
//! \param num_nodes The number of nodes in the input graph.
inline size_t ThetaPrime(ssize_t x, double epsilonPrime, double l, size_t k,
                         size_t num_nodes, mpi_omp_parallel_tag &&) {
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  return (ThetaPrime(x, epsilonPrime, l, k, num_nodes, omp_parallel_tag{}) /
          world_size) +
         1;
}

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

template <typename GraphTy, typename ConfTy, typename RRRGeneratorTy,
          typename diff_model_tag, typename execution_tag>
std::pair<std::vector<unsigned int>, int> MartigaleRound(
  TimerAggregator &timeAggregator,
  ssize_t thetaPrime,  
  TransposeRRRSets<GraphTy>& tRRRSets,
  std::vector<RRRset<GraphTy>> RR,
  RRRsetAllocator<typename GraphTy::vertex_type> allocator,
  const GraphTy &G, 
  const ConfTy &CFG, 
  RRRGeneratorTy &generator, 
  IMMExecutionRecord &record,
  diff_model_tag &&model_tag, 
  execution_tag &&ex_tag, 
  std::vector<int> vertexToProcess,
  int world_size, 
  int world_rank,
  CommunicationEngine<GraphTy> cEngine
) 
{
  double f;
  // figure out how to exacly generate this number (remove rounding error)
  ssize_t localThetaPrime = (thetaPrime / world_size) + 1;

  std::pair<std::vector<unsigned int>, int> globalSeeds;

  size_t delta = localThetaPrime - RR.size();
  record.ThetaPrimeDeltas.push_back(delta);

  auto timeRRRSets = measure<>::exec_time([&]() {
    RR.insert(RR.end(), delta, RRRset<GraphTy>(allocator));
    auto begin = RR.end() - delta;

    // std::cout << "generating transpose RRR sets within the martigale loop" << std::endl;
    // within this function, there is a logical bug with counting RRRset IDs
    timeAggregator.samplingTimer.startTimer();
    GenerateTransposeRRRSets(tRRRSets, G, generator, begin, RR.end(), record,
                    std::forward<diff_model_tag>(model_tag),
                    std::forward<execution_tag>(ex_tag));
    timeAggregator.samplingTimer.endTimer();

    // std::cout << "start linearization rank = " << world_rank << std::endl;
    LinearizedSetsSize* setSize = cEngine.count(tRRRSets, vertexToProcess, world_size);
    int* data = cEngine.linearize(
      tRRRSets, 
      vertexToProcess, 
      /// TODO: Most likely a memroy leak, free this data at some point
      *(cEngine.buildPrefixSum(setSize->countPerProcess.data(), world_size)), 
      setSize->count, 
      world_size
    );

    // std::cout << "linearizing data, rank = " << world_rank << std::endl;
    std::unordered_map<int, std::unordered_set<int>>* aggregateSets = new std::unordered_map<int, std::unordered_set<int>>();

    timeAggregator.allToAllTimer.startTimer();
    // auto start = std::chrono::high_resolution_clock::now();
    cEngine.getProcessSpecificVertexRRRSets(*aggregateSets, data, setSize->countPerProcess.data(), world_size, localThetaPrime);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << " ------- time for single all to all = " << (end - start).count() << " ------- " << std::endl;
    timeAggregator.allToAllTimer.endTimer();

    delete setSize;

    // std::cout << "first max_cover_greedy, rank = " << world_rank << std::endl;
    timeAggregator.max_k_localTimer.startTimer();

    /// TODO: This needs to output a std::pair<unordered_set<unsigned int>, int> instead of a vector. This will allow you to do constant time 
    ///   lookup in the next stage when you linearize these results. 
    MaxKCoverEngine localKCoverEngine((int)CFG.k);
    std::pair<std::vector<unsigned int>, int> localSeeds = localKCoverEngine.useLazyGreedy(*aggregateSets)->run_max_k_cover(*aggregateSets, thetaPrime*2);
    std::unordered_set<unsigned int> localSeedsSet(localSeeds.first.begin(), localSeeds.first.end());

    // end = std::chrono::high_resolution_clock::now();
    timeAggregator.max_k_localTimer.endTimer();

    // linearize the winning k seeds from each local process (vertexID: {RRR set IDs}) and send them to process 1 (reduce operation). 
    // std::cout << "linearize, rank = " << world_rank << std::endl;

    /// TODO: Linear Local Seeds is generated from the output of maxKCoverEngine.max_cover + the aggregateSets. Currently you are sending all 
    ///   of your data to the global node (node 0) instead of just the top k seeds from each vertex. This is why it's so slow. 
    std::pair<int, int*> linearLocalSeeds = cEngine.linearizeLocalSeeds(*aggregateSets, localSeedsSet);

    // std::cout << "global aggregation, rank = " << world_rank << std::endl;
    // mpi_gather for rank 1. Gathers the size of all linearLocalSeeds to send
    timeAggregator.allGatherTimer.startTimer();
    std::pair<int, int*> globalAggregation = cEngine.aggregateAggregateSets(linearLocalSeeds.first, world_size, linearLocalSeeds.second);
    timeAggregator.allGatherTimer.endTimer();

    int* aggregatedSeeds = globalAggregation.second;
    int totalData = globalAggregation.first;

    /// TODO: This will become more expensive with each martigale iteration (the rest of a martigale iteration 
    // does not become more expensive with each iteration due to the progressive data generation). This means that
    // when k is very large this can become very expensive. Look into ways to optimize the aggregate key seleciton.
    if (world_rank == 0) {
      // delinearize the m * k vertex: unordered_set data. 
      std::unordered_map<int, std::unordered_set<int>> bestKMSeeds;

      // std::cout << "getting all local seeds seeds, rank = " << world_rank << std::endl;
      cEngine.aggregateLocalKSeeds(bestKMSeeds, aggregatedSeeds, totalData);

      // run max-k-cover on the aggregated data for k seeds
      // std::cout << "calculating global seeds, rank = " << world_rank << std::endl;
      timeAggregator.max_k_globalTimer.startTimer();
      MaxKCoverEngine globalKCoverEngine((int)CFG.k);
      globalSeeds = globalKCoverEngine.useStochasticGreedy((double)CFG.epsilon)->useLazyGreedy(bestKMSeeds)->run_max_k_cover(bestKMSeeds, thetaPrime * 2);
      // end = std::chrono::high_resolution_clock::now();
      timeAggregator.max_k_globalTimer.endTimer();
    }    

    // free statements to prevent memory leaks.
    delete aggregatedSeeds;
    delete linearLocalSeeds.second;
    delete aggregateSets;
  });
  
  record.ThetaEstimationGenerateRRR.push_back(timeRRRSets);

  auto timeMostInfluential = measure<>::exec_time([&]() { });

  record.ThetaEstimationMostInfluential.push_back(timeMostInfluential);

  // std::cout << "returning global seeds" << std::endl;
  return globalSeeds;
}


//! Collect a set of Random Reverse Reachable set.
//!
//! \tparam GraphTy The type of the input graph.
//! \tparam RRRGeneratorTy The type of the RRR generator.
//! \tparam diff_model_tag Type-Tag to selecte the diffusion model.
//! \tparam execution_tag Type-Tag to select the execution policy.
//!
//! \param G The input graph.  The graph is transoposed.
//! \param k The size of the seed set.
//! \param epsilon The parameter controlling the approximation guarantee.
//! \param l Parameter usually set to 1.
//! \param generator The rrr sets generator.
//! \param record Data structure storing timing and event counts.
//! \param model_tag The diffusion model tag.
//! \param ex_tag The execution policy tag.
template <typename GraphTy, typename ConfTy, typename RRRGeneratorTy,
          typename diff_model_tag, typename execution_tag>
std::pair<std::vector<unsigned int>, int> TransposeSampling(
  TransposeRRRSets<GraphTy>& tRRRSets,
  const GraphTy &G, const ConfTy &CFG, double l,
  RRRGeneratorTy &generator, IMMExecutionRecord &record,
  diff_model_tag &&model_tag, execution_tag &&ex_tag, 
  std::vector<int> vertexToProcess,
  int world_size, int world_rank
) {
  using vertex_type = typename GraphTy::vertex_type;
  size_t k = CFG.k;
  double epsilon = CFG.epsilon;
  
  // sqrt(2) * epsilon
  double epsilonPrime = 1.4142135623730951 * epsilon;

  // each iteration of the martigale loop does
    // sampling
    // communication 
    // max-k-cover 

  // if martigale success return best local seeds

  double LB = 0;
  #if defined ENABLE_MEMKIND
  RRRsetAllocator<vertex_type> allocator(libmemkind::kinds::DAX_KMEM_PREFERRED);
  #elif defined ENABLE_METALL
  RRRsetAllocator<vertex_type> allocator =  metall_manager_instance().get_allocator();
  #else
  RRRsetAllocator<vertex_type> allocator;
  #endif
  std::vector<RRRset<GraphTy>> RR;

  auto start = std::chrono::high_resolution_clock::now();
  size_t thetaPrime = 0;

  TimerAggregator timeAggregator;

  CommunicationEngine<GraphTy> cEngine;

  // martingale loop
  for (ssize_t x = 1; x < std::log2(G.num_nodes()); ++x) {
    // Equation 9
    // thetaPrime == number of RRRsets globally
    ssize_t thetaPrime = ThetaPrime(x, epsilonPrime, l, k, G.num_nodes(),
                                    std::forward<execution_tag>(ex_tag));

    // if f(s) meets threashold return k seeds, all processes return null except for rank 
    // std::cout << "martigale round with theta = " << (int)thetaPrime << std::endl;
    std::pair<std::vector<unsigned int>, int> seeds = MartigaleRound(
      timeAggregator,
      thetaPrime, tRRRSets, RR, allocator, G, 
      CFG, generator, record,
      std::forward<diff_model_tag>(model_tag),
      std::forward<execution_tag>(ex_tag),
      vertexToProcess, world_size, world_rank,
      cEngine
    );

    // f is the fraction of RRRsets covered by the seeds / the total number of RRRSets (in the current iteration of the martigale loop)
    // this has to be a global value, if one process succeeds and another fails it will get stuck in communication (the algorithm will fail). 
    double f;

    if (world_rank == 0) {
      f = double(seeds.second) / thetaPrime;
    }
    // mpi_broadcast f(s)

    timeAggregator.broadcastTimer.startTimer();
    cEngine.distributeF(&f);
    timeAggregator.broadcastTimer.endTimer();

    // std::cout << "seeds.second: (covered RRRSet IDs) = " << seeds.second << " , thetaPrme: " << thetaPrime << " , f = " << f << std::endl;
    if (f >= std::pow(2, -x)) {
      // std::cout << "Fraction " << f << std::endl;
      LB = (G.num_nodes() * f) / (1 + epsilonPrime);
      spdlog::get("console")->info("Lower Bound {}", LB);
      break;
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  size_t theta = Theta(epsilon, l, k, LB, G.num_nodes());
  
  if (theta == 0) {
    theta = thetaPrime;
  }

  size_t localTheta = (theta / world_size) + 1;

  record.ThetaEstimationTotal = end - start;

  record.Theta = theta;
  spdlog::get("console")->info("Theta {}", theta);

  record.GenerateRRRSets = measure<>::exec_time([&]() {
    if (localTheta > RR.size()) {
      size_t final_delta = localTheta - RR.size();
      RR.insert(RR.end(), final_delta, RRRset<GraphTy>(allocator));

      auto begin = RR.end() - final_delta;
      timeAggregator.samplingTimer.startTimer();
      GenerateTransposeRRRSets(tRRRSets, G, generator, begin, RR.end(), record,
                      std::forward<diff_model_tag>(model_tag),
                      std::forward<execution_tag>(ex_tag));
      timeAggregator.samplingTimer.endTimer();
    }
  });

  // std::cout << "getting best seeds globally with theta = " << (int)theta << std::endl;
  std::pair<std::vector<unsigned int>, int> bestSeeds = MartigaleRound(
    timeAggregator,
    theta, tRRRSets, RR, allocator, G, 
    CFG, generator, record,
    std::forward<diff_model_tag>(model_tag),
    std::forward<execution_tag>(ex_tag),
    vertexToProcess, world_size, world_rank,
    cEngine
  );

  // std::cout << "total communication time: " << cEngine.getCommunicationTime() << std::endl;
  // std::cout << "total sample time: " << totalSampleTime << std::endl;
  // std::cout << "total max_k_cover time: " << totalCoverTime << std::endl;

  std::cout << "Samping time: " << timeAggregator.samplingTimer.resolveTimer() << std::endl;
  std::cout << "Max-k-cover local process time: " << timeAggregator.max_k_localTimer.resolveTimer() << std::endl;
  std::cout << "Max-k-cover global process time: " << timeAggregator.max_k_globalTimer.resolveTimer() << std::endl;
  std::cout << "AllGather time: " << timeAggregator.allGatherTimer.resolveTimer() << std::endl;
  std::cout << "AlltoAll time: " << timeAggregator.allToAllTimer.resolveTimer() << std::endl;
  std::cout << "Broadcast time: " << timeAggregator.broadcastTimer.resolveTimer() << std::endl;

  return bestSeeds;
}


//! The IMM algroithm for Influence Maximization (MPI specialization).
//!
//! \tparam GraphTy The type of the input graph.
//! \tparam PRNG The type of the parallel random number generator.
//! \tparam diff_model_tag Type-Tag to selecte the diffusion model.
//!
//! \param G The input graph.  The graph is transoposed.
//! \param k The size of the seed set.
//! \param epsilon The parameter controlling the approximation guarantee.
//! \param l Parameter usually set to 1.
//! \param gen The parallel random number generator.
//! \param model_tag The diffusion model tag.
//! \param ex_tag The execution policy tag.
template <typename GraphTy, typename ConfTy, typename diff_model_tag,
          typename GeneratorTy, typename ExTagTrait>
auto GREEDI(const GraphTy &G, const ConfTy &CFG, double l, GeneratorTy &gen,
         IMMExecutionRecord &record, 
         diff_model_tag &&model_tag, ExTagTrait &&) {
  using vertex_type = typename GraphTy::vertex_type;

  l = l * (1 + 1 / std::log2(G.num_nodes()));

  // get mapping of which local process handles each vertex
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // TOOD: randomize the bucketing of vertices to processes, use a coinflip algorithm
  //  that chooses a process for each vertex as a loop linearly scans through the 
  //  verticesPerProcess vector.
  std::vector<int> vertexToProcess(G.num_nodes(), -1);
  int verticesPerProcess = G.num_nodes() / world_size + 1;
  for (int i = 0; i < G.num_nodes(); i++) {
    vertexToProcess[i] = i % world_size;
  }

  TransposeRRRSets<GraphTy>* tRRRSets = new TransposeRRRSets<GraphTy>(G.num_nodes());
  
  // get local tRRRSets
  // std::cout << "entering transposeSampling" << std::endl;
  std::pair<std::vector<unsigned int>, int> seeds = mpi::TransposeSampling(
    *tRRRSets,
    G, CFG, l, gen, record,
    std::forward<diff_model_tag>(model_tag),
    typename ExTagTrait::generate_ex_tag{},
    vertexToProcess,
    world_size, world_rank
  );

  delete tRRRSets;

  return seeds.first;
}

}  // namespace mpi
}  // namespace ripples

#endif  // RIPPLES_MPI_IMM_H
