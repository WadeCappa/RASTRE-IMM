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

#ifndef RIPPLES_GRAPH_H
#define RIPPLES_GRAPH_H

#include <algorithm>
#include <cstddef>
#include <map>
#include <numeric>
#include <vector>

#include <ripples/utility.h>

namespace ripples {

//! \brief Forward Direction Graph loading policy.
//!
//! \tparam VertexTy The type of the vertex in the graph.
template <typename VertexTy>
struct ForwardDirection {
  //! \brief Edge Source
  //!
  //! \tparam ItrTy Edge iterator type.
  //!
  //! \param itr Iterator to the current edge.
  //! \return The source of the egde to be loaded in the graph.
  template <typename ItrTy, typename MapTy>
  static VertexTy Source(ItrTy itr, const MapTy &m) {
    return m.find(itr->source)->second;
  }

  //! \brief Edge Destination
  //!
  //! \tparam ItrTy Edge iterator type.
  //!
  //! \param itr Iterator to the current edge.
  //! \return The destination of the egde to be loaded in the graph.
  template <typename ItrTy, typename MapTy>
  static VertexTy Destination(ItrTy itr, const MapTy & m) {
    return m.find(itr->destination)->second;
  }
};

//! \brief Backward Direction Graph loading policy.
//!
//! \tparam VertexTy The type of the vertex in the graph.
template <typename VertexTy>
struct BackwardDirection {
  //! \brief Edge Source
  //!
  //! \tparam ItrTy Edge iterator type.
  //!
  //! \param itr Iterator to the current edge.
  //! \return The source of the egde to be loaded in the graph.
  template <typename ItrTy, typename MapTy>
  static VertexTy Source(ItrTy itr, const MapTy & m) {
    return m.find(itr->destination)->second;
  }

  //! \brief Edge Destination
  //!
  //! \tparam ItrTy Edge iterator type.
  //!
  //! \param itr Iterator to the current edge.
  //! \return The destination of the egde to be loaded in the graph.
  template <typename ItrTy, typename MapTy>
  static VertexTy Destination(ItrTy itr, const MapTy & m) {
    return m.find(itr->source)->second;
  }
};

//! \brief A weighted edge.
//!
//! \tparam VertexTy The integer type representing a vertex of the graph.
//! \tparam WeightTy The type representing the weight on the edge.
template <typename VertexTy, typename WeightTy>
struct Edge {
  //! The integer type representing vertices in the graph.
  using vertex_type = VertexTy;
  //! The type representing weights on the edges of the graph.
  using weight_type = WeightTy;
  //! The source of the edge.
  VertexTy source;
  //! The destination of the edge.
  VertexTy destination;
  //! The weight on the edge.
  WeightTy weight;

  bool operator==(const Edge & O) const {
    return O.source == this->source
        && O.destination == this->destination
        && O.weight == this->weight;
  }
};

//! \brief The Graph data structure.
//!
//! A graph in CSR format.  The construction method takes care of projecting the
//! vercites in the input edge list into a contiguous space [0; N[ of integers
//! in order to build the CSR representation.  However, the data structure
//! stores a map that allows to project back the IDs into the original space.
//!
//! \tparam VertexTy The integer type representing a vertex of the graph.
//! \tparam WeightTy The type representing the weight on the edge.
template <typename VertexTy, typename WeightTy,
          typename DirectionPolicy = ForwardDirection<VertexTy>>
class Graph {
 public:
  //! The size type.
  using size_type = size_t;
  //! The type of an edge in the graph.
  using edge_type = Edge<VertexTy, WeightTy>;
  //! The integer type representing vertices in the graph.
  using vertex_type = VertexTy;
  //! The type representing the weights on the edges of the graph.
  using edge_weight_type = WeightTy;

  //! \brief The edges stored in the CSR.
  struct DestinationTy {
    VertexTy vertex;  //!< The destination of an edges.
    WeightTy weight;  //!< The edge weight.

    bool operator==(const DestinationTy& O) const {
      return this->vertex == O.vertex && this->weight == O.weight;
    }
  };

  //! \brief The neighborhood of a vertex.
  class Neighborhood {
   public:
    //! Construct the neighborhood.
    //!
    //! \param B The begin of the neighbor list.
    //! \param E The end of the neighbor list.
    Neighborhood(DestinationTy *B, DestinationTy *E) : begin_(B), end_(E) {}

    //! Begin of the neighborhood.
    //! \return an iterator to the begin of the neighborhood.
    DestinationTy *begin() const { return begin_; }
    //! End of the neighborhood.
    //! \return an iterator to the begin of the neighborhood.
    DestinationTy *end() const { return end_; }

   private:
    DestinationTy *begin_;
    DestinationTy *end_;
  };

  //! Empty Graph Constructor.
  Graph()
      : numNodes(0),
        numEdges(0),
        index(nullptr),
        edges(nullptr),
        idMap(),
        reverseMap() {}

  Graph(const Graph &O)
      : numNodes(O.numNodes)
      , numEdges(O.numEdges)
      , idMap(O.idMap)
      , reverseMap(O.reverseMap) {
    edges = new DestinationTy[numEdges];
    index = new DestinationTy * [numNodes + 1];
    std::copy(O.edges, O.edges + numEdges, edges);
    std::transform(O.index, O.index + numNodes + 1, index,
                   [=](const DestinationTy * p) -> DestinationTy * {
                     return edges + p - O.index;
                   });
  }

  Graph & operator=(const Graph &O) {
    numNodes = O.numNodes;
    numEdges = O.numEdges;
    idMap = O.idMap;
    reverseMap = O.reverseMap;

    edges = new DestinationTy[numEdges];
    index = new DestinationTy * [numNodes + 1];
    std::copy(O.edges, O.edges + numEdges, edges);
    std::transform(O.index, O.index + numNodes + 1, index,
                   [=](const DestinationTy * p) -> DestinationTy * {
                     return edges + p - O.index;
                   });
  }

  //! Move constructor.
  //! \param O The graph to be moved.
  Graph(Graph &&O) {
    std::swap(numNodes, O.numNodes);
    std::swap(numEdges, O.numEdges);
    std::swap(index, O.index);
    std::swap(edges, O.edges);
    idMap = std::move(O.idMap);
    reverseMap = std::move(O.reverseMap);
  }

  //! Move assignment operator.
  //! \param O The graph to be moved.
  //! \return a reference to the destination graph.
  Graph &operator=(Graph &&O) {
    std::swap(numNodes, O.numNodes);
    std::swap(numEdges, O.numEdges);
    std::swap(index, O.index);
    std::swap(edges, O.edges);
    idMap = std::move(O.idMap);
    reverseMap = std::move(O.reverseMap);
    return *this;
  }

  //! Reload from binary constructor.
  //!
  //! \tparam FStream The type of the input stream.
  //!
  //! \param FS The binary stream containing the graph dump.
  template <typename FStream>
  Graph(FStream &FS) {
    load_binary(FS);
  }

  //! \brief Constructor.
  //!
  //! Build a Graph from a sequence of edges.  The vertex identifiers are
  //! projected over the integer interval [0;N[.  The data structure stores
  //! conversion maps to move fro the internal representation of the vertex IDs
  //! to the original input representation.
  //!
  //! \tparam EdgeIterator The iterator type used to visit the input edge list.
  //!
  //! \param begin The start of the edge list.
  //! \param end The end of the edge list.
  template <typename EdgeIterator>
  Graph(EdgeIterator begin, EdgeIterator end) {
    for (auto itr = begin; itr != end; ++itr) {
      idMap[itr->source];
      idMap[itr->destination];
    }

    size_t num_nodes = idMap.size();
    size_t num_edges = std::distance(begin, end);

    index = new DestinationTy *[num_nodes + 1];
    edges = new DestinationTy[num_edges];

#pragma omp parallel for simd
    for (size_t i = 0; i < num_nodes + 1; ++i) {
      index[i] = edges;
    }

#pragma omp parallel for simd
    for (size_t i = 0; i < num_edges; ++i) {
      edges[i] = DestinationTy();
    }

    numNodes = num_nodes;
    numEdges = num_edges;

    VertexTy currentID{0};
    reverseMap.resize(numNodes);
    for (auto itr = std::begin(idMap), end = std::end(idMap); itr != end;
         ++itr) {
      reverseMap[currentID] = itr->first;
      itr->second = currentID;
      currentID++;
    }

    for (auto itr = begin; itr != end; ++itr) {
      index[DirectionPolicy::Source(itr, idMap) + 1] += 1;
    }

    for (size_t i = 1; i <= num_nodes; ++i) {
      index[i] += index[i - 1] - edges;
    }

    std::vector<DestinationTy *> ptrEdge(index, index + num_nodes);
    for (auto itr = begin; itr != end; ++itr) {
      *ptrEdge[DirectionPolicy::Source(itr, idMap)] = {
        DirectionPolicy::Destination(itr, idMap), itr->weight};
      ++ptrEdge[DirectionPolicy::Source(itr, idMap)];
    }
  }

  //! \brief Destuctor.
  ~Graph() {
    delete[] index;
    delete[] edges;
  }

  //! Returns the out-degree of a vertex.
  //! \param v The input vertex.
  //! \return the in-degree of vertex v in input.
  size_t degree(VertexTy v) const { return index[v + 1] - index[v]; }

  //! Returns the neighborhood of a vertex.
  //! \param v The input vertex.
  //! \return  a range containing the out-neighbors of the vertex v in input.
  Neighborhood neighbors(VertexTy v) const {
    return Neighborhood(index[v], index[v + 1]);
  }

  //! The number of nodes in the Graph.
  //! \return The number of nodes in the Graph.
  size_t num_nodes() const { return numNodes; }

  //! The number of edges in the Graph.
  //! \return The number of edges in the Graph.
  size_t num_edges() const { return numEdges; }

  //! Convert a list of vertices from the interal representation to the original
  //! input representation.
  //!
  //! \tparam Itr The iterator type of the input sequence of vertex IDs.
  //! \tparam OutputItr The iterator type of the output sequence.
  //!
  //! \param b The begin of the input vertex IDs sequence.
  //! \param e The end of the input vertex IDs sequence.
  //! \param o The start of the output vertex IDs sequence.
  template <typename Itr, typename OutputItr>
  void convertID(Itr b, Itr e, OutputItr o) const {
    using value_type = typename Itr::value_type;
    std::transform(b, e, o, [&](const value_type &v) -> value_type {
        return reverseMap.at(v);
    });
  }

  //! Convert a vertex from the interal representation to the original input
  //! representation.
  //!
  //! \param v The input vertex ID.
  //! \return The original vertex ID in the input representation.
  vertex_type convertID(const vertex_type v) const {
    return reverseMap.at(v);
  }

  //! Convert a list of vertices from the original input edge list
  //! representation to the internal vertex representation.
  //!
  //! \tparam Itr The iterator type of the input sequence of vertex IDs.
  //! \tparam OutputItr The iterator type of the output sequence.
  //!
  //! \param b The begin of the input vertex IDs sequence.
  //! \param e The end of the input vertex IDs sequence.
  //! \param o The start of the output vertex IDs sequence.
  template <typename Itr, typename OutputItr>
  void transformID(Itr b, Itr e, OutputItr o) const {
    using value_type = typename Itr::value_type;
    std::transform(b, e, o, [this](const value_type &v) -> value_type {
        return transformID(v);
    });
  }

  vertex_type transformID(const vertex_type v) const {
      auto itr = idMap.find(v);
      if (itr != idMap.end())
        return itr->second;
      else
        throw "Bad node";
  }

  //! Dump the internal representation to a binary stream.
  //!
  //! \tparam FStream The type of the output stream
  //!
  //! \param FS The ouput file stream.
  template <typename FStream>
  void dump_binary(FStream &FS) const {
    uint64_t num_nodes = htole64(numNodes);
    uint64_t num_edges = htole64(numEdges);
    FS.write(reinterpret_cast<const char *>(&num_nodes), sizeof(uint64_t));
    FS.write(reinterpret_cast<const char *>(&num_edges), sizeof(uint64_t));

    sequence_of<VertexTy>::dump(FS, reverseMap.begin(), reverseMap.end());

    using relative_index =
        typename std::iterator_traits<DestinationTy *>::difference_type;
    std::vector<relative_index> relIndex(numNodes + 1, 0);
    std::transform(index, index + numNodes + 1, relIndex.begin(),
                   [=](DestinationTy *v) -> relative_index {
                     return std::distance(edges, v);
                   });
    sequence_of<relative_index>::dump(FS, relIndex.begin(), relIndex.end());
    sequence_of<DestinationTy>::dump(FS, edges, edges + numEdges);
  }

 private:
  static constexpr bool isForward =
      std::is_same<DirectionPolicy, ForwardDirection<VertexTy>>::value;
  using transposed_direction =
      typename std::conditional<isForward, BackwardDirection<VertexTy>,
                                ForwardDirection<VertexTy>>::type;
  using transposed_type =
      Graph<vertex_type, edge_weight_type, transposed_direction>;

  friend transposed_type;
 public:

  //! Get the transposed graph.
  //! \return the transposed graph.
  transposed_type get_transpose() const {
    using out_dest_type = typename transposed_type::DestinationTy;
    transposed_type G;
    G.numEdges = numEdges;
    G.numNodes = numNodes;
    G.reverseMap = reverseMap;
    G.idMap = idMap;
    G.index = new out_dest_type *[numNodes + 1];
    G.edges = new out_dest_type[numEdges];

    std::fill(G.index, G.index + numNodes + 1, nullptr);

    std::for_each(edges, edges + numEdges,
                  [&](const DestinationTy &d) { ++G.index[d.vertex + 1]; });

    G.index[0] = G.edges;
    std::partial_sum(G.index, G.index + numNodes + 1, G.index,
                     [](out_dest_type *a, out_dest_type *b) -> out_dest_type * {
                       size_t sum = reinterpret_cast<size_t>(a) +
                                    reinterpret_cast<size_t>(b);
                       return reinterpret_cast<out_dest_type *>(sum);
                     });

    std::vector<out_dest_type *> destPointers(G.index, G.index + numNodes);

    for (vertex_type v = 0; v < numNodes; ++v) {
      for (auto u : neighbors(v)) {
        *destPointers[u.vertex] = {v, u.weight};
        destPointers[u.vertex]++;
      }
    }

    return G;
  }

 private:
  template <typename FStream>
  void load_binary(FStream &FS) {
    if (!FS.is_open()) throw "Bad things happened!!!";

    FS.read(reinterpret_cast<char *>(&numNodes), sizeof(numNodes));
    FS.read(reinterpret_cast<char *>(&numEdges), sizeof(numEdges));

    numNodes = le64toh(numNodes);
    numEdges = le64toh(numEdges);

    reverseMap.resize(numNodes);
    FS.read(reinterpret_cast<char *>(reverseMap.data()),
            reverseMap.size() * sizeof(VertexTy));

    sequence_of<VertexTy>::load(reverseMap.begin(), reverseMap.end(),
                                reverseMap.begin());

    for (VertexTy i = 0; i < numNodes; ++i) idMap[reverseMap[i]] = i;

    index = new DestinationTy *[numNodes + 1];
    edges = new DestinationTy[numEdges];

    FS.read(reinterpret_cast<char *>(index),
            (numNodes + 1) * sizeof(ptrdiff_t));

    sequence_of<DestinationTy *>::load(index, index + numNodes + 1, index);

    std::transform(index, index + numNodes + 1, index,
                   [=](DestinationTy *v) -> DestinationTy * {
                     return reinterpret_cast<ptrdiff_t>(v) + edges;
                   });

    FS.read(reinterpret_cast<char *>(edges), numEdges * sizeof(DestinationTy));
    sequence_of<DestinationTy>::load(edges, edges + numEdges, edges);
  }

  DestinationTy **index;
  DestinationTy *edges;

  std::map<VertexTy, VertexTy> idMap;
  std::vector<VertexTy> reverseMap;

  size_t numNodes;
  size_t numEdges;
};


template <typename BwdGraphTy, typename FwdGraphTy>
auto getCommunitiesSubgraphs(const FwdGraphTy & Gf,
                             const std::vector<typename FwdGraphTy::vertex_type> & communityVector) {
  using vertex_type = typename FwdGraphTy::vertex_type;
  size_t num_communities = *std::max_element(communityVector.begin(), communityVector.end()) + 1;
  std::vector<BwdGraphTy> communities(num_communities);

  using EdgeTy = typename FwdGraphTy::edge_type;
  std::vector<std::vector<EdgeTy>> edge_lists(num_communities);
  for (typename FwdGraphTy::vertex_type src = 0; src < Gf.num_nodes(); ++src) {
    vertex_type original_src = Gf.convertID(src);

    vertex_type community_src = communityVector[original_src - 1];
    for (auto e : Gf.neighbors(src)) {
      vertex_type original_dst = Gf.convertID(e.vertex);

      vertex_type community_dst = communityVector[original_dst - 1];
      if (community_dst == community_src) {
        edge_lists[community_src].push_back({original_src, original_dst, e.weight});
      }
    }
  }

  for (size_t i = 0; i < num_communities; ++i) {
    communities[i] = BwdGraphTy(edge_lists[i].begin(), edge_lists[i].end());
  }

  return communities;
}

}  // namespace ripples

#endif  // RIPPLES_GRAPH_H
