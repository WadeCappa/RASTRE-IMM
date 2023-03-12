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

#include "catch2/catch.hpp"

#include "ripples/generate_rrr_sets.h"
#include "ripples/graph.h"
#include "ripples/imm.h"
#include "ripples/mpi/communication_engine.h"
#include "ripples/max_k_cover.h"

#include <unordered_set>
#include <unordered_map>
#include <iostream>

#include "trng/lcg64.hpp"

using EdgeT = ripples::Edge<uint32_t, float>;
int received_data[] = { 
    0, 25, 84, 32, 3, 0, -1, 
    1, 91, 86, 62, 57, 55, 13, 42, 9, 67, 37, 66, 92, 33, 4, 6, 64, 97, 71, 12, 41, 1, 30, 16, 45, 78, 19, 48, 20, 80, 21, 50, 46, 23, 83, 24, 53, 29, 54, -1, 
    2, 92, 79, 75, 73, 91, 62, 86, 57, 55, 54, 29, 78, 49, 23, 1, 12, 97, 64, 6, 19, 67, 8, 37, 45, 16, 13, 42, 18, -1, 
    3, 86, 73, 62, 55, 83, 54, 29, 1, 64, 97, 68, 16, 42, 67, 9, 22, 8, 37, 50, 20, -1, 
    4, 78, 57, 62, 49, 42, 37, 86, 73, 8, 67, -1, 
    0, 80, 70, 53, 66, -1, 
    1, 56, 30, 48, 80, 41, -1, 
    2, 85, 80, 30, 66, -1, 
    3, 57, 37, 8, 67, 6, -1, 
    4, 44, 55, 29, 16, 23, 10, 9, 86, 21, 99, 65, -1, 
    0, 86, 26, 42, 2, -1, 
    1, 53, -1, 
    2, -1, 
    3, 94, -1, 
    4, 86, 78, 62, 20, 42, 69, -1, 
};

int receivedDataSizes[3] = { 111, 39, 22 };
int expectedOutputSizesPerVertex[5] = { 13, 44, 32, 26, 27 };
int p = 3; 
int RRRIDsPerProcess = 100;

int k = 2; 
int theta = p * RRRIDsPerProcess;


SCENARIO("Alltoall and alltoall_v have already been run, the local process must aggregate the data.", "[aggregate_tRRRsets]") {
  GIVEN("Serialized data received from alltoall_v") {

    using destination_type = ripples::WeightedDestination<uint32_t, float>;
    using GraphBwd = ripples::Graph<uint32_t, destination_type,
                                    ripples::BackwardDirection<uint32_t>>;

    WHEN("alltoall_v and alltoall have been successfully run") {
        CommunicationEngine<GraphBwd> cEngine;
        std::unordered_map<int, std::unordered_set<int>> local_tRRRsets;
        cEngine.aggregateTRRRSets(local_tRRRsets, received_data, receivedDataSizes, p, RRRIDsPerProcess);

        THEN("Total data count remains the same, but data has been aggregated") {
            REQUIRE(local_tRRRsets.size() == 5);

            for (auto i : local_tRRRsets)
            {
                std::cout << "vertexID: " << i.first << std::endl;
                for (auto j : i.second)
                {
                    std::cout << j << " ";
                }
                std::cout << std::endl;
                std::cout << "expected size: " << expectedOutputSizesPerVertex[i.first] << " real size: " << i.second.size() << std::endl;
                REQUIRE(expectedOutputSizesPerVertex[i.first] == i.second.size());
            }
        }

        THEN("max k-cover is run locally") {
            MaxKCoverEngine kCoverEngine;
            std::vector<int> res = kCoverEngine.max_cover_lazy_greedy(local_tRRRsets, k, theta);
            REQUIRE(res.size() == k);
            for (const auto & vertex : res) {
                std::cout << vertex << std::endl;
            }
        }
    }
  }
}
