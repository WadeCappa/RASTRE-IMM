#include <vector>
#include <set>
#include <unordered_set>
#include <mutex>
#include <iostream>

template <typename GraphTy>
class TransposeRRRSets
{
    public:
    std::vector<std::pair<std::mutex, std::set<typename GraphTy::vertex_type>>> sets;

    public:
    TransposeRRRSets(int numberOfVertices) 
        : sets(numberOfVertices)
    {
    }

    ~TransposeRRRSets() 
    {
    }

    void addToSet(int index, typename GraphTy::vertex_type vertex) 
    {
        sets[index].first.lock(); 
        sets[index].second.insert(vertex);
        sets[index].first.unlock();
    }

    auto getBegin() 
    {
        return sets.begin();
    }

    auto getEnd() 
    {
        return sets.end();
    }
};

