//===------------------------------------------------------------*- C++ -*-===//
//
// Copyright 2018 Battelle Memorial Institute
//
//===----------------------------------------------------------------------===//

#ifndef IM_FIND_MOST_INFLUENTIAL_H
#define IM_FIND_MOST_INFLUENTIAL_H

#include <algorithm>
#include <queue>
#include <vector>

#include <omp.h>

#include "im/utility.h"

namespace im {

template <typename ItrTy1, typename ItrTy2>
ItrTy2 swap_ranges(ItrTy1 B, ItrTy1 E, ItrTy2 O, sequential_tag) {
  return std::swap_ranges(B, E, O);
}

template <typename ItrTy1, typename ItrTy2>
ItrTy2 swap_ranges(ItrTy1 B, ItrTy1 E, ItrTy2 O, omp_parallel_tag) {
  size_t toBeSwaped = std::distance(B, E);
#pragma omp parallel for
  for (size_t i = 0; i < toBeSwaped; ++i) {
    std::iter_swap(B + i, O + i);
  }
  return O + toBeSwaped;
}

template <typename ItrTy, typename UnaryPredicate>
ItrTy partition(ItrTy B, ItrTy E, UnaryPredicate P, sequential_tag) {
  return std::partition(B, E, P);
}

namespace {

template <typename ItrTy, typename ex_tag = omp_parallel_tag>
struct PartitionIndices {
  ItrTy begin;
  ItrTy end;
  ItrTy pivot;

  PartitionIndices() : begin{0}, end{0}, pivot{0} {}

  bool operator==(const PartitionIndices &O) const {
    return this->begin == O.begin && this->end == O.end &&
           this->pivot == O.pivot;
  }

  PartitionIndices &operator+=(const PartitionIndices &O) {
    if (*this == PartitionIndices()) {
      this->begin = O.begin;
      this->end = O.end;
      this->pivot = O.pivot;
      return *this;
    }

    if (O == PartitionIndices()) {
      return *this;
    }

    if (this->pivot == this->begin && O.pivot == O.begin) {
      this->end = O.end;
      return *this;
    } else if (this->pivot == this->end && O.pivot == O.end) {
      this->end = O.end;
      this->pivot = O.pivot;
      return *this;
    }

    if (std::distance(this->pivot, this->end) <
        std::distance(O.begin, O.pivot)) {
      size_t toBeMoved = std::distance(this->pivot, this->end);
      swap_ranges(this->pivot, this->end, std::prev(O.pivot, toBeMoved),
                  ex_tag{});
      this->pivot = std::prev(O.pivot, toBeMoved);
    } else {
      this->pivot = std::swap_ranges(O.begin, O.pivot, this->pivot);
    }
    this->end = O.end;

    return *this;
  }
};

}  // namespace

template <typename ItrTy, typename UnaryPredicate>
ItrTy partition(ItrTy B, ItrTy E, UnaryPredicate P, omp_parallel_tag) {
  size_t num_threads(1);

#pragma omp single
  num_threads = omp_get_max_threads();

  std::vector<PartitionIndices<ItrTy>> indices(num_threads);

#pragma omp parallel
  {
    size_t num_elements = std::distance(B, E);
    size_t threadnum = omp_get_thread_num(), numthreads = omp_get_num_threads();
    size_t low = num_elements * threadnum / numthreads,
           high = num_elements * (threadnum + 1) / numthreads;

    indices[threadnum].begin = B + low;
    indices[threadnum].end = std::min(E, B + high);
    indices[threadnum].pivot =
        std::partition(indices[threadnum].begin, indices[threadnum].end, P);
  }

#pragma omp parallel
  {
#pragma omp single
    {
      for (size_t j = 1; j < num_threads; j <<= 1) {
        for (size_t i = 0; (i + j) < num_threads; i += j * 2) {
#pragma omp task
          { indices[i] += indices[i + j]; }
        }
#pragma omp taskwait
      }
    }
  }

  return indices[0].pivot;
}

//! \brief Count the occurrencies of vertices in the RRR sets.
//!
//! \tparam InItr The input sequence iterator type.
//! \tparam OutItr The output sequence iterator type.
//!
//! \param in_begin The begin of the sequence of RRR sets.
//! \param in_end The end of the sequence of RRR sets.
//! \param out_begin The begin of the sequence storing the counters for each
//! vertex.
//! \param out_end The end of the sequence storing the counters for each vertex.
template <typename InItr, typename OutItr>
void CountOccurrencies(InItr in_begin, InItr in_end, OutItr out_begin,
                       OutItr out_end, sequential_tag &&) {
  using rrr_set_type = typename std::iterator_traits<InItr>::value_type;
  using vertex_type = typename rrr_set_type::value_type;
  for (; in_begin != in_end; ++in_begin) {
    std::for_each(in_begin->begin(), in_begin->end(),
                  [&](const vertex_type v) { *(out_begin + v) += 1; });
  }
}

//! \brief Count the occurrencies of vertices in the RRR sets.
//!
//! \tparam InItr The input sequence iterator type.
//! \tparam OutItr The output sequence iterator type.
//!
//! \param in_begin The begin of the sequence of RRR sets.
//! \param in_end The end of the sequence of RRR sets.
//! \param out_begin The begin of the sequence storing the counters for each
//! vertex.
//! \param out_end The end of the sequence storing the counters for each vertex.
template <typename InItr, typename OutItr>
void CountOccurrencies(InItr in_begin, InItr in_end, OutItr out_begin,
                       OutItr out_end, omp_parallel_tag &&) {
  using rrr_set_type = typename std::iterator_traits<InItr>::value_type;
  using vertex_type = typename rrr_set_type::value_type;

#pragma omp parallel
  {
    size_t num_elements = std::distance(out_begin, out_end);
    size_t threadnum = omp_get_thread_num(), numthreads = omp_get_num_threads();
    vertex_type low = num_elements * threadnum / numthreads,
                high = num_elements * (threadnum + 1) / numthreads;

    for (auto itr = in_begin; itr != in_end; ++itr) {
      auto begin = std::lower_bound(itr->begin(), itr->end(), low);
      auto end = std::upper_bound(begin, itr->end(), high - 1);
      std::for_each(begin, end,
                    [&](const vertex_type v) { *(out_begin + v) += 1; });
    }
  }
}

//! \brief Initialize the Heap storage.
//!
//! \tparam InItr The input sequence iterator type.
//! \tparam OutItr The output sequence iterator type.
//!
//! \param in_begin The begin of the sequence of vertex counters.
//! \param in_end The end of the sequence of vertex counters.
//! \param out_begin The begin of the sequence used as storage in the Heap.
//! \param out_end The end of the sequence used as storage in the Heap.
template <typename InItr, typename OutItr>
void InitHeapStorage(InItr in_begin, InItr in_end, OutItr out_begin,
                     OutItr out_end, sequential_tag &&) {
  using value_type = typename std::iterator_traits<OutItr>::value_type;
  using vertex_type = typename value_type::first_type;

  for (vertex_type v = 0; in_begin != in_end; ++in_begin, ++v, ++out_begin) {
    *out_begin = {v, *in_begin};
  }
}

//! \brief Initialize the Heap storage.
//!
//! \tparam InItr The input sequence iterator type.
//! \tparam OutItr The output sequence iterator type.
//!
//! \param in_begin The begin of the sequence of vertex counters.
//! \param in_end The end of the sequence of vertex counters.
//! \param out_begin The begin of the sequence used as storage in the Heap.
//! \param out_end The end of the sequence used as storage in the Heap.
template <typename InItr, typename OutItr>
void InitHeapStorage(InItr in_begin, InItr in_end, OutItr out_begin,
                     OutItr out_end, omp_parallel_tag &&) {
  using value_type = typename std::iterator_traits<OutItr>::value_type;
  using vertex_type = typename value_type::first_type;

#pragma omp parallel for
  for (vertex_type v = 0; v < std::distance(in_begin, in_end); ++v) {
    *(out_begin + v) = {v, *(in_begin + v)};
  }
}

//! \brief Update the coverage counters.
//!
//! \tparam VertexTy The type of the vertices.
//! \tparam RRRsetsTy The type storing RRR sets.
//! \tparam RemovedVectorTy The type of the vector storing removed RRR sets.
//! \tparam VertexCoverageVectorTy The type of the vector storing counters.
//!
//! \param v The chosen vertex.
//! \param RRRsets The sequence of RRRsets.
//! \param removed The vector storing covered/uncovered flags for the RRRsets.
//! \param vertexCoverage The vector storing the counters to be updated.
template <typename VertexTy, typename RRRsetsTy, typename RemovedVectorTy,
          typename VertexCoverageVectorTy>
void UpdateCounters(const VertexTy v, const RRRsetsTy &RRRsets,
                    RemovedVectorTy &removed,
                    VertexCoverageVectorTy &vertexCoverage, sequential_tag &&) {
  for (size_t i = 0; i < RRRsets.size(); ++i) {
    if (removed[i]) continue;

    if (std::binary_search(RRRsets[i].begin(), RRRsets[i].end(), v)) {
      removed[i] = true;
      for (size_t j = 0; j < RRRsets[i].size(); ++j) {
        vertexCoverage[RRRsets[i][j]] -= 1;
      }
    }
  }
}

//! \brief Update the coverage counters.
//!
//! \tparam VertexTy The type of the vertices.
//! \tparam RRRsetsTy The type storing RRR sets.
//! \tparam RemovedVectorTy The type of the vector storing removed RRR sets.
//! \tparam VertexCoverageVectorTy The type of the vector storing counters.
//!
//! \param v The chosen vertex.
//! \param RRRsets The sequence of RRRsets.
//! \param removed The vector storing covered/uncovered flags for the RRRsets.
//! \param vertexCoverage The vector storing the counters to be updated.
template <typename VertexTy, typename RRRsetsTy, typename RemovedVectorTy,
          typename VertexCoverageVectorTy>
void UpdateCounters(const VertexTy v, const RRRsetsTy &RRRsets,
                    RemovedVectorTy &removed,
                    VertexCoverageVectorTy &vertexCoverage,
                    omp_parallel_tag &&) {
  for (size_t i = 0; i < RRRsets.size(); ++i) {
    if (removed[i]) continue;

    if (std::binary_search(RRRsets[i].begin(), RRRsets[i].end(), v)) {
      removed[i] = true;

#pragma omp parallel for
      for (size_t j = 0; j < RRRsets[i].size(); ++j) {
        vertexCoverage[RRRsets[i][j]] -= 1;
      }
    }
  }
}

//! \brief Select k seeds starting from the a list of Random Reverse
//! Reachability Sets.
//!
//! \tparam GraphTy The graph type.
//! \tparam RRRset The type storing Random Reverse Reachability Sets.
//! \tparam execution_tag The execution policy.
//!
//! \param G The input graph.
//! \param k The size of the seed set.
//! \param RRRsets A vector of Random Reverse Reachability sets.
//! \param ex_tag The execution policy tag.
//!
//! \return a pair where the size_t is the number of RRRset covered and
//! the set of vertices selected as seeds.
template <typename GraphTy, typename RRRset, typename execution_tag>
auto FindMostInfluentialSet(const GraphTy &G, size_t k,
                            const std::vector<RRRset> &RRRsets,
                            execution_tag &&ex_tag) {
  using vertex_type = typename GraphTy::vertex_type;

  std::vector<uint32_t> vertexCoverage(G.num_nodes(), 0);

  auto cmp = [](std::pair<vertex_type, size_t> &a,
                std::pair<vertex_type, size_t> &b) {
    return a.second < b.second;
  };
  using priorityQueue =
      std::priority_queue<std::pair<vertex_type, size_t>,
                          std::vector<std::pair<vertex_type, size_t>>,
                          decltype(cmp)>;

  std::vector<std::pair<vertex_type, size_t>> queue_storage(G.num_nodes());

  CountOccurrencies(RRRsets.begin(), RRRsets.end(), vertexCoverage.begin(),
                    vertexCoverage.end(), std::forward<execution_tag>(ex_tag));

  InitHeapStorage(vertexCoverage.begin(), vertexCoverage.end(),
                  queue_storage.begin(), queue_storage.end(),
                  std::forward<execution_tag>(ex_tag));

  priorityQueue queue(cmp, std::move(queue_storage));

  std::vector<typename GraphTy::vertex_type> result;
  result.reserve(k);

  std::vector<char> removed(RRRsets.size(), false);
  size_t uncovered = RRRsets.size();

  while (result.size() < k && uncovered != 0) {
    auto element = queue.top();
    queue.pop();

    if (element.second > vertexCoverage[element.first]) {
      element.second = vertexCoverage[element.first];
      queue.push(element);
      continue;
    }

    uncovered -= element.second;

    UpdateCounters(element.first, RRRsets, removed, vertexCoverage,
                   std::forward<execution_tag>(ex_tag));

    result.push_back(element.first);
  }

  double f = double(RRRsets.size() - uncovered) / RRRsets.size();

  return std::make_pair(f, result);
}

}  // namespace im

#endif  // IM_FIND_MOST_INFLUENTIAL_H
