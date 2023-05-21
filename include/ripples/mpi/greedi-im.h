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

#include "ripples/mpi/streaming_engine.h"

namespace ripples {
namespace mpi {

// seeds = MartigaleRound(
//   timeAggregator,
//   thetaPrime, tRRRSets, RR, allocator, G, 
//   CFG, gen, record,
//   std::forward<diff_model_tag>(model_tag),
//   std::forward<execution_tag>(ex_tag),
//   vertexToProcess, world_size, world_rank,
//   cEngine,
//   aggregateSets,
//   old_sampling_sizes
// );

template <
  typename GraphTy, 
  typename ConfTy, 
  typename RRRGeneratorTy,
  typename diff_model_tag, 
  typename execution_tag
>
class RanDIMM
{
  private:
  const GraphTy &G;
  const ConfTy &CFG;
  ripples::omp_parallel_tag &ex_tag;
  diff_model_tag &model_tag;
  RRRGeneratorTy &gen;
  IMMExecutionRecord &record;
  const CommunicationEngine<GraphTy> &cEngine;

  const double l;
  size_t thetaPrime = 0;

  std::map<int, std::vector<int>> aggregateSets;
  std::vector<size_t> old_sampling_sizes;
  size_t RR_sets;
  RRRsetAllocator<typename GraphTy::vertex_type> allocator;
  TransposeRRRSets<GraphTy> tRRRSets;
  const std::vector<int> vertexToProcess;

  TimerAggregator &timeAggregator;

  RRRsetAllocator<typename GraphTy::vertex_type> GetAllocator()
  {
    #if defined ENABLE_MEMKIND
      return allocator(libmemkind::kinds::DAX_KMEM_PREFERRED);
    #elif defined ENABLE_METALL
      return metall_manager_instance().get_allocator();
    #else
      return RRRsetAllocator<typename GraphTy::vertex_type>();
    #endif
  }

std::pair<std::vector<unsigned int>, int> MartigaleRound(
  ssize_t thetaPrime
) 
{
  const size_t localThetaPrime = (thetaPrime / this->cEngine.GetSize()) + 1;

  std::pair<std::vector<unsigned int>, size_t> approximated_solution;

  auto timeRRRSets = measure<>::exec_time([&]() {

    size_t delta = localThetaPrime - this->RR_sets;
    record.ThetaPrimeDeltas.push_back(delta);

    spdlog::get("console")->info("sampling...");

    timeAggregator.samplingTimer.startTimer();
    
    GenerateTransposeRRRSets(
      tRRRSets, this->RR_sets, delta, G, this->gen, record,
      std::forward<diff_model_tag>(model_tag),
      std::forward<execution_tag>(ex_tag)
    );

    this->RR_sets += delta;

    timeAggregator.samplingTimer.endTimer();    

    spdlog::get("console")->info("AlltoAll...");

    timeAggregator.allToAllTimer.startTimer();  

    cEngine.AggregateThroughAllToAll(
      this->tRRRSets,
      this->vertexToProcess,
      this->old_sampling_sizes,
      localThetaPrime,
      this->aggregateSets
    );
  
    timeAggregator.allToAllTimer.endTimer();

    spdlog::get("console")->info("seed selection...");

    int kprime = int(CFG.alpha * (double)CFG.k);

    if (CFG.use_streaming == true) 
    {
      if (this->cEngine.GetRank() == 0) 
      {
        StreamingRandGreedIEngine streamingEngine((int)CFG.k, kprime, thetaPrime*2, (double)CFG.epsilon_2, this->cEngine.GetSize() - 1);
        approximated_solution = streamingEngine.Stream(&timeAggregator);
      }
      else 
      {
        timeAggregator.totalSendTimer.startTimer();
        
        StreamingMaxKCover<GraphTy> localKCoverEngine((int)CFG.k, kprime, timeAggregator, this->cEngine);
        localKCoverEngine.useLazyGreedy(this->aggregateSets);
        localKCoverEngine.run_max_k_cover(this->aggregateSets, thetaPrime*2);

        timeAggregator.totalSendTimer.endTimer();
      }
    }
    else 
    {
      timeAggregator.max_k_localTimer.startTimer();

      std::pair<std::vector<unsigned int>, ssize_t> localSeeds;

      {
        MaxKCover<GraphTy> localKCoverEngine((int)CFG.k, kprime, timeAggregator);
        localSeeds = localKCoverEngine.useLazyGreedy(this->aggregateSets).run_max_k_cover(this->aggregateSets, thetaPrime*2);
      }

      timeAggregator.max_k_localTimer.endTimer();

      spdlog::get("console")->info("all gather...");
      timeAggregator.allGatherTimer.startTimer();

      std::vector<unsigned int> globalAggregation;
      size_t totalData = this->cEngine.GatherPartialSolutions(localSeeds, this->aggregateSets, globalAggregation);

      timeAggregator.allGatherTimer.endTimer();

      if (this->cEngine.GetRank() == 0) {
        std::map<int, std::vector<int>> bestKMSeeds;

        spdlog::get("console")->info("unpacking local seeds in global process...");

        timeAggregator.allGatherTimer.startTimer();
        std::vector<std::pair<unsigned int, std::vector<unsigned int>*>>* local_seeds = cEngine.aggregateLocalKSeeds(bestKMSeeds, globalAggregation.data(), totalData);

        timeAggregator.allGatherTimer.endTimer();

        spdlog::get("console")->info("global max_k_cover...");
        timeAggregator.max_k_globalTimer.startTimer();

        {
          MaxKCover<GraphTy> globalKCoverEngine((int)CFG.k, kprime, timeAggregator);
          approximated_solution = globalKCoverEngine.useLazyGreedy(bestKMSeeds).run_max_k_cover(bestKMSeeds, thetaPrime * 2);
        }

        // end = std::chrono::high_resolution_clock::now();

        for (const auto & s: *local_seeds)
        {
          if (s.first > approximated_solution.second)
          {
            approximated_solution.second = s.first;
            approximated_solution.first = *(s.second);
          }
        }

        timeAggregator.max_k_globalTimer.endTimer();

        for (auto & s : *local_seeds)
        {
          delete s.second;
        }

        delete local_seeds;
      }
    }    

  });
  
  record.ThetaEstimationGenerateRRR.push_back(timeRRRSets);
  auto timeMostInfluential = measure<>::exec_time([&]() { });
  record.ThetaEstimationMostInfluential.push_back(timeMostInfluential);

  std::cout << "rank " << this->cEngine.GetRank() << " returning global seeds" << std::endl;
  return approximated_solution;
}

  public:
  RanDIMM
  (
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
    this->allocator = this->GetAllocator();
    this->RR_sets = 0;

    for (size_t i = 0; i < vertexToProcess.size(); i++)
    {
      if (vertexToProcess[i] == this->cEngine.GetRank())
      {
        this->aggregateSets.insert({i, std::vector<int>()});
      }
    }
  }

  auto ApproximateSeedSet(const size_t theta)
  {
    return this->MartigaleRound(
      theta
    );
  }

  void OutputDiagnosticData()
  {
    if (this->CFG.use_streaming == true)
    {
      std::cout << " --- SHARED --- " << std::endl; 
      std::cout << "Samping time: " << this->timeAggregator.samplingTimer.resolveTimer() << std::endl;
      std::cout << "AlltoAll time: " << this->timeAggregator.allToAllTimer.resolveTimer() << std::endl;
      std::cout << "Receive Broadcast: " << this->timeAggregator.broadcastTimer.resolveTimer() << std::endl;

      std::cout << " --- SENDER --- " << std::endl; 
      std::cout << "Select Next Seed: " << this->timeAggregator.max_k_localTimer.resolveTimer() << std::endl;
      std::cout << "Send Next Seed: " << this->timeAggregator.sendTimer.resolveTimer() << std::endl;
      std::cout << "Total Send Time: " << this->timeAggregator.totalSendTimer.resolveTimer() << std::endl;
      
      std::cout << " --- RECEIVER --- " << std::endl; 
      std::cout << "Initialize Buckets: " << this->timeAggregator.initBucketTimer.resolveTimer() << std::endl;
      std::cout << "Receive Next Seed: " << this->timeAggregator.receiveTimer.resolveTimer() << std::endl;
      std::cout << "Insert Into Buckets: " << this->timeAggregator.max_k_globalTimer.resolveTimer() << std::endl;
      std::cout << "Handling received data (inserting into matrix and copying from buffer): " << this->timeAggregator.processingReceiveTimer.resolveTimer() << std::endl; 
      std::cout << "Atomic Update (receiver side): " << this->timeAggregator.atomicUpdateTimer.resolveTimer() << std::endl; 
      std::cout << "Total Global Streaming Time: " << this->timeAggregator.totalGlobalStreamTimer.resolveTimer() << std::endl;
    } 
    else  
    {
      std::cout << " --- SHARED --- " << std::endl; 
      std::cout << "Samping time: " << this->timeAggregator.samplingTimer.resolveTimer() << std::endl;
      std::cout << "f score Broadcast time: " << this->timeAggregator.broadcastTimer.resolveTimer() << std::endl;
      std::cout << "AlltoAll time: " << this->timeAggregator.allToAllTimer.resolveTimer() << std::endl;
      std::cout << "AllGather time: " << this->timeAggregator.allGatherTimer.resolveTimer() << std::endl;

      std::cout << " --- LOCAL --- " << std::endl; 
      std::cout << "Local max-cover time: " << this->timeAggregator.max_k_localTimer.resolveTimer() << std::endl;

      std::cout << " --- GLOBAL --- " << std::endl; 
      std::cout << "Global max-cover time: " << this->timeAggregator.max_k_globalTimer.resolveTimer() << std::endl;
    }
  }
};

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
// template <typename GraphTy, typename ConfTy, typename RRRGeneratorTy,
//           typename diff_model_tag, typename execution_tag>
// std::pair<std::vector<unsigned int>, int> TransposeSampling(
//   TransposeRRRSets<GraphTy>& tRRRSets,
//   const GraphTy &G, const ConfTy &CFG, double l,
//   RRRGeneratorTy &generator, IMMExecutionRecord &record,
//   diff_model_tag &&model_tag, execution_tag &&ex_tag, 
//   std::vector<int>& vertexToProcess,
//   int world_size, int world_rank
// ) {
//   using vertex_type = typename GraphTy::vertex_type;

// }

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
          typename RRRGeneratorTy, typename ExTagTrait>
auto GREEDI(
  const GraphTy &G, 
  const ConfTy &CFG, 
  double l_value, 
  RRRGeneratorTy &gen,
  IMMExecutionRecord &record, 
  diff_model_tag &&model_tag, 
  ExTagTrait &&
) 
{
  using vertex_type = typename GraphTy::vertex_type;

  // no sequential version available
  using execution_tag = ripples::omp_parallel_tag;

  ////// INITIALIZATION //////

  CommunicationEngine<GraphTy> cEngine = CommunicationEngineBuilder<GraphTy>::BuildCommunicationEngine();
  execution_tag ex_tag = typename ExTagTrait::generate_ex_tag{};

  TimerAggregator timeAggregator;

  // std::cout << "before initialization" << std::endl;
  // std::cout << "number of nodes before functions " << G.num_nodes() << std::endl;

  RanDIMM<GraphTy, ConfTy, RRRGeneratorTy, diff_model_tag, execution_tag> randimm(
    G, CFG, ex_tag, model_tag, l_value, gen, record, cEngine, timeAggregator
  );

  double LB = 0;
  double epsilonPrime = 1.4142135623730951 * CFG.epsilon;

  size_t thetaPrime = 0;

  std::pair<std::vector<unsigned int>, int> seeds;
  
  std::map<int, std::vector<int>> aggregateSets;
  
  ////// MARTINGALE LOOP //////

  std::cout << "starting martingale loop" << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  for (ssize_t x = 1; x < std::log2(G.num_nodes()); ++x) 
  {
    // Equation 9
    ssize_t thetaPrime = ThetaPrime(
      x, epsilonPrime, l_value, CFG.k, G.num_nodes(), 
      std::forward<execution_tag>(ex_tag)
    );

    seeds = randimm.ApproximateSeedSet(thetaPrime);

    std::cout << "finished iteration " << x << " and aquired utility of " << seeds.second << std::endl;

    // f is the fraction of RRRsets covered by the seeds / the total number of RRRSets (in the current iteration of the martigale loop)
    // this has to be a global value, if one process succeeds and another fails it will get stuck in communication (the algorithm will fail). 
    double f;

    std::cout << "thetaprime: " << thetaPrime << std::endl;

    if (cEngine.GetRank() == 0) {
      f = (double)(seeds.second) / thetaPrime;
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

  size_t theta = Theta(CFG.epsilon, l_value, CFG.k, LB, G.num_nodes());
  
  if (theta == 0) {
    theta = thetaPrime;
  }

  record.ThetaEstimationTotal = end - start;

  record.Theta = theta;
  spdlog::get("console")->info("Theta {}", theta);

  std::pair<std::vector<vertex_type>, int> bestSeeds;

  if (theta <= thetaPrime) 
  {
    bestSeeds = seeds; 
  }
  else 
  {
    bestSeeds = randimm.ApproximateSeedSet(theta);
  }

  randimm.OutputDiagnosticData();

  return seeds.first;
}

}  // namespace mpi
}  // namespace ripples

#endif  // RIPPLES_MPI_IMM_H
