#include <vector>
#include <set>
#include <unordered_set>
#include <mutex>
#include <iostream>
#include "ripples/bitmask.h"

template <typename GraphTy>
class TransposeRRRSets
{
    public:
    std::vector<std::pair<std::mutex, std::vector<size_t>>> sets;

    public:
    TransposeRRRSets(size_t numberOfVertices) 
        : sets(numberOfVertices)
    {
    }

    ~TransposeRRRSets() 
    {
    }

    void addToSet(int vertex, size_t RRRId) 
    {
        sets[vertex].first.lock(); 
        sets[vertex].second.push_back(RRRId);
        sets[vertex].first.unlock();
    }

    unsigned int calculateInfluence(const std::vector<unsigned int> &seed_set, const size_t theta) const {
        ripples::Bitmask<int> covered(theta);

        #pragma omp parallel for
        for (const auto & vertex : seed_set) {
            for (const auto & RRR_set : this->sets[vertex].second) {
                covered.set(RRR_set);
            }
        }

        return covered.popcount();
    }

    void GetSetSizes(std::vector<unsigned int>& setSizes)
    {
        for (int i = 0; i < setSizes.size() && i < this->sets.size(); i++)
        {
            setSizes[i] = (unsigned(this->sets[i].second.size()));
        }
    }

    void RemoveDuplicates()
    {
        #pragma omp parallel for
        for (int i = 0; i < this->sets.size(); i++)
        {
            int sizeBefore = this->sets[i].second.size();

            std::unordered_set<size_t> seen;
            for (const auto & r : this->sets[i].second)
            {
                seen.insert(r);
            }

            this->sets[i].second.assign(seen.begin(), seen.end());

            if (sizeBefore != this->sets[i].second.size())
            {
                std::cout << "size difference; " << sizeBefore -  this->sets[i].second.size() << std::endl;
            }
        }
    }
};

