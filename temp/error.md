In file included from ../../include/ripples/generate_rrr_sets.h:55,
                 from ../../test/rrr_set_generation.cc:45:
../../include/ripples/transposeRRRSets.h: In instantiation of ‘void TransposeRRRSets<GraphTy>::addToSet(int, typename GraphTy::vertex_type) [with GraphTy = ripples::Graph<unsigned int, ripples::WeightedDestination<unsigned int, float>, ripples::BackwardDirection<unsigned int> >; typename GraphTy::vertex_type = unsigned int]’:
../../include/ripples/generate_rrr_sets.h:190:20:   required from ‘void ripples::AddTransposeRRRSet(TransposeRRRSets<GraphTy>&, const GraphTy&, typename GraphTy::vertex_type, PRNGeneratorTy&, RRRset<GraphTy>&, diff_model_tag&&) [with GraphTy = Graph<unsigned int, WeightedDestination<unsigned int, float>, BackwardDirection<unsigned int> >; PRNGeneratorTy = trng::lcg64; diff_model_tag = independent_cascade_tag; typename GraphTy::vertex_type = unsigned int; RRRset<GraphTy> = std::vector<unsigned int>]’
../../include/ripples/streaming_rrr_generator.h:198:25:   required from ‘void ripples::CPUWalkWorker<GraphTy, PRNGeneratorTy, ItrTy, diff_model_tag>::batchTranspose(TransposeRRRSets<GraphTy>&, ItrTy, ItrTy) [with GraphTy = ripples::Graph<unsigned int, ripples::WeightedDestination<unsigned int, float>, ripples::BackwardDirection<unsigned int> >; PRNGeneratorTy = trng::lcg64; ItrTy = __gnu_cxx::__normal_iterator<std::vector<unsigned int>*, std::vector<std::vector<unsigned int> > >; diff_model_tag = ripples::independent_cascade_tag]’
../../include/ripples/streaming_rrr_generator.h:179:7:   required from ‘void ripples::CPUWalkWorker<GraphTy, PRNGeneratorTy, ItrTy, diff_model_tag>::svc_transpose_loop(std::atomic<long unsigned int>&, TransposeRRRSets<GraphTy>&, ItrTy, ItrTy) [with GraphTy = ripples::Graph<unsigned int, ripples::WeightedDestination<unsigned int, float>, ripples::BackwardDirection<unsigned int> >; PRNGeneratorTy = trng::lcg64; ItrTy = __gnu_cxx::__normal_iterator<std::vector<unsigned int>*, std::vector<std::vector<unsigned int> > >; diff_model_tag = ripples::independent_cascade_tag]’
../../include/ripples/streaming_rrr_generator.h:170:8:   required from here
../../include/ripples/transposeRRRSets.h:28:16: error: no match for ‘operator+’ (operand types are ‘std::vector<std::pair<std::mutex, std::unordered_set<unsigned int, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<unsigned int> > >, std::allocator<std::pair<std::mutex, std::unordered_set<unsigned int, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<unsigned int> > > > >’ and ‘int’)
   28 |         (*sets + index).first.lock();
      |         ~~~~~~~^~~~~~~~
In file included from /usr/include/c++/12/string:47,
                 from ../../../../.conan/data/catch2/2.13.7/_/_/package/5ab84d6acfe1f23c4fae0ab88f26e3a396351ac9/include/catch2/catch.hpp:475,
                 from ../../test/rrr_set_generation.cc:43:
/usr/include/c++/12/bits/stl_iterator.h:630:5: note: candidate: ‘template<class _Iterator> std::reverse_iterator<_Iterator> std::operator+(typename reverse_iterator<_Iterator>::difference_type, const reverse_iterator<_Iterator>&)’
  630 |     operator+(typename reverse_iterator<_Iterator>::difference_type __n,
      |     ^~~~~~~~
/usr/include/c++/12/bits/stl_iterator.h:630:5: note:   template argument deduction/substitution failed:
../../include/ripples/transposeRRRSets.h:28:16: note:   mismatched types ‘const std::reverse_iterator<_Iterator>’ and ‘int’
   28 |         (*sets + index).first.lock();
      |         ~~~~~~~^~~~~~~~
/usr/include/c++/12/bits/stl_iterator.h:1786:5: note: candidate: ‘template<class _Iterator> std::move_iterator<_IteratorL> std::operator+(typename move_iterator<_IteratorL>::difference_type, const move_iterator<_IteratorL>&)’
 1786 |     operator+(typename move_iterator<_Iterator>::difference_type __n,
      |     ^~~~~~~~
/usr/include/c++/12/bits/stl_iterator.h:1786:5: note:   template argument deduction/substitution failed:
../../include/ripples/transposeRRRSets.h:28:16: note:   mismatched types ‘const std::move_iterator<_IteratorL>’ and ‘int’
