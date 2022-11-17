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
#include "ripples/communication_engine.h"

#include <unordered_set>
#include <iostream>

#include "trng/lcg64.hpp"

using EdgeT = ripples::Edge<uint32_t, float>;
std::vector<int> data_process_1 = { 
    0, 25, 84, 32, 3, 0, -1, 
    1, 91, 86, 62, 57, 55, 13, 42, 9, 67, 37, 66, 92, 33, 4, 6, 64, 97, 71, 12, 41, 1, 30, 16, 45, 78, 19, 48, 20, 80, 21, 50, 46, 23, 83, 24, 53, 29, 54, -1, 
    2, 92, 79, 75, 73, 91, 62, 86, 57, 55, 54, 29, 78, 49, 23, 1, 12, 97, 64, 6, 19, 67, 8, 37, 45, 16, 13, 42, 18, -1, 
    3, 86, 73, 62, 55, 83, 54, 29, 1, 64, 97, 68, 16, 42, 67, 9, 22, 8, 37, 50, 20, -1, 
    4, 78, 57, 62, 49, 42, 37, 86, 73, 8, 67, -1, 
    5, 80, 70, 53, 66, -1, 
    6, 56, 30, 48, 80, 41, -1, 
    7, 85, 80, 30, 66, -1, 
    8, 57, 37, 8, 67, 6, -1, 
    9, 44, 55, 29, 16, 23, 10, 9, 86, 21, 99, 65, -1, 
};

std::vector<int> data_process_2 = { 
    0, 86, 26, 42, 2, -1, 
    1, 53, -1, 
    2, -1, 
    3, 94, -1, 
    4, 86, 78, 62, 20, 42, 69, -1, 
    5, 88, 86, 16, -1, 
    6, 61, 86, 21, 42, 51, 38, -1, 
    7, 80, 30, -1, 
    8, 75, 52, 5, -1, 
    9, 82, 59, 27, 81, 42, 10, 9, 86, 34, 58, 71, -1, 
};

std::vector<int> vertex_mapping = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};


SCENARIO("Alltoall and alltoallV between two processes.", "[transpose_rrrsets]") {
  GIVEN("Serialized data and a mapping for each vertex to a process.") {


    WHEN("alltoall and alltoallv are used AFTER local processes generate tRRRSets.") {

    }
  }
}
