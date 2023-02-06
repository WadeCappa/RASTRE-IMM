#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <mutex>
#include <iostream>
#include "bitmask.h"
#include <utility>
#include <queue>
#include <cstdlib>
#include <bits/stdc++.h>
#include <cmath>

class MaxKCoverEngine 
{
private:

    class NextMostInfluentialFinder
    { 
    protected:
        std::vector<unsigned int>* vertex_subset;
        std::unordered_map<int, std::unordered_set<int>*>* allSets;
        size_t subset_size;

    public:
        virtual ~NextMostInfluentialFinder()
        {
            // std::cout << "deallocating base class..." << std::endl;
            for (const auto & l : *(this->allSets))
            {
                delete l.second;
            }
            delete this->allSets;
        }

        virtual int findNextInfluential (
            std::vector<unsigned int>& seedSet,
            int current_K_index,
            ripples::Bitmask<int>& covered,
            int theta
        ) = 0;

        virtual NextMostInfluentialFinder* setSubset (
            std::vector<unsigned int>* subset_of_selection_sets,
            size_t subset_size
        ) = 0;

        virtual NextMostInfluentialFinder* reloadSubset () = 0;
    };

    class LazyGreedy : public NextMostInfluentialFinder
    {
    private:
        template <typename idTy>
        struct CompareMaxHeap {

            bool operator()(const std::pair<idTy, std::unordered_set<idTy>*> a,
                            const std::pair<idTy, std::unordered_set<idTy>*> b) {
                return a.second->size() < b.second->size();
            }
        };

        CompareMaxHeap<int> cmp;            
        std::vector<std::pair<int, std::unordered_set<int>*>>* heap;

    public:
        LazyGreedy(std::unordered_map<int, std::unordered_set<int>>& data)
        {
            this->allSets = new std::unordered_map<int, std::unordered_set<int>*>();

            for (const auto & l : data)
            {
                this->allSets->insert({ l.first, new std::unordered_set<int>(l.second.begin(), l.second.end()) });
            }
        }

        ~LazyGreedy()
        {
            // std::cout << "deallocating Lazy-Greedy finder ..." << std::endl;
            delete heap;
        }

        NextMostInfluentialFinder* setSubset(std::vector<unsigned int>* subset_of_selection_sets, size_t subset_size) override
        {
            this->subset_size = subset_size;
            this->heap = new std::vector<std::pair<int, std::unordered_set<int>*>>(this->subset_size);

            generateQueue(subset_of_selection_sets, subset_size);
            this->vertex_subset = subset_of_selection_sets;
            return this;
        }

        NextMostInfluentialFinder* reloadSubset () override 
        {
            generateQueue(this->vertex_subset, this->subset_size);
            return this;
        }

        int findNextInfluential(
            std::vector<unsigned int>& seedSet,
            int current_K_index,
            ripples::Bitmask<int>& covered,
            int theta
        ) override
        {
            int totalCovered;
            std::pair<int, std::unordered_set<int>*> l = this->heap->front();
            std::pop_heap(this->heap->begin(), this->heap->end(), this->cmp);

            std::unordered_set<int> temp;

            // remove RR IDs from l that are already covered. 
            for (int e: *(l.second)) {
                if (e > theta || e < 0) {
                    std::cout << "e is greater than theta, e = " << e << " , theta = " << theta << std::endl;
                }
                else if (covered.get(e)) {
                    temp.insert(e);
                }
            }            

            for (const int e : temp) 
            {
                l.second->erase(e); 
            }
            
            // calculate marginal gain
            auto marginal_gain = l.second->size();

            // calculate utiluty of next best
            auto r = this->heap->front();
            
            // if marginal of l is better than r's utility, l is the current best     
            if (marginal_gain >= r.second->size()) {
                seedSet[current_K_index] = l.first;
                
                for (int e : *(l.second)) {
                    if (e > theta || e < 0) {
                        std::cout << "e is greater than theta, e = " << e << " , theta = " << theta << std::endl;
                    }
                    else if (!covered.get(e)) {
                        totalCovered++;
                        covered.set(e);
                    }
                }
            }
            // else push l's marginal into the heap 
            else {
                this->heap->push_back(l);
                std::push_heap(this->heap->begin(), this->heap->end(), this->cmp);
                return findNextInfluential(
                    seedSet, current_K_index, covered, theta
                );
            }

            return totalCovered;
        }

        private:
        void generateQueue(std::vector<unsigned int>* subset_of_selection_sets, size_t subset_size)
        {
            for (int i = 0; i < subset_size; i++)
            {
                this->heap->at(i) = std::make_pair(subset_of_selection_sets->at(i), this->allSets->at(subset_of_selection_sets->at(i)));
            }

            std::make_heap(this->heap->begin(), this->heap->end(), this->cmp);
        }
    };

    class NaiveGreedy : public NextMostInfluentialFinder
    {
    public:
        NaiveGreedy(std::unordered_map<int, std::unordered_set<int>>& data) 
        {
            this->allSets = new std::unordered_map<int, std::unordered_set<int>*>();

            for (const auto & l : data)
            {
                this->allSets->insert({ l.first, new std::unordered_set<int>(l.second.begin(), l.second.end()) });
            }
        }

        ~NaiveGreedy()
        {
            // std::cout << "deallocating Naive-Greedy finder ..." << std::endl;
        }

        NextMostInfluentialFinder* setSubset(std::vector<unsigned int>* subset_of_selection_sets, size_t subset_size) override
        {
            this->vertex_subset = subset_of_selection_sets;
            this->subset_size = subset_size;
            return this;
        } 

        NextMostInfluentialFinder* reloadSubset () override 
        {
            return this;
        }

        int findNextInfluential(
            std::vector<unsigned int>& seedSet,
            int current_K_index,
            ripples::Bitmask<int>& covered,
            int theta
        ) override
        {
            int max = 0;
            int max_key = -1;
            int totalCovered = 0;

            for ( int i = 0; i < this->subset_size; i++ )
            {
                int vertex = this->vertex_subset->at(i);
                if (this->allSets->find(vertex) != this->allSets->end() && this->allSets->at(vertex)->size() > max)
                {
                    max = this->allSets->at(vertex)->size();
                    max_key = vertex;
                }
            }
            seedSet[current_K_index] = max_key;

            for (int e: *(this->allSets->at(max_key))) {
                if (!covered.get(e)) {
                    covered.set(e);
                    totalCovered++;
                }
            }

            #pragma omp parallel 
            {
                # pragma omp for schedule(static)
                for( int i = 0; i < this->subset_size; i++ ) {
                    if (this->allSets->find(this->vertex_subset->at(i)) != this->allSets->end()) 
                    {
                        auto RRRSets = this->allSets->at(this->vertex_subset->at(i));

                        std::set<int> temp;
                        if (this->vertex_subset->at(i) != max_key) {
                            for (int e: *RRRSets) {
                                if (covered.get(e)) {
                                    temp.insert(e);
                                }
                            }
                            for (int e: temp) {
                                this->allSets->at(this->vertex_subset->at(i))->erase(e); 
                            }

                        }
                    }
                }
            }

            this->allSets->erase(max_key);
            return totalCovered;
        }
    };

    int k;
    double epsilon;
    bool usingStochastic = false;
    NextMostInfluentialFinder* finder = 0;

    void reorganizeVertexSet(std::vector<unsigned int>* vertices, size_t size, std::vector<unsigned int> seedSet)
    {
        // for i from 0 to n−2 do
        //     j ← random integer such that i ≤ j < n
        //     exchange a[i] and a[j]
        std::unordered_set<int> seeds(seedSet.begin(), seedSet.end());

        for (int i = 0; i < size; i++) 
        {
            int j = getRandom(i, vertices->size()-1);
            while (seeds.find(vertices->at(j)) != seeds.end())
            {
                j = getRandom(i, vertices->size()-1);
            }
            std::swap(vertices->at(i), vertices->at(j));
        }
    }

    int getRandom(int min, int max)
    {
        return min + (rand() % static_cast<int>(max - min + 1));
    }

    size_t getSubsetSize(size_t n, int k, double epsilon)
    {
        // new set R = (n/k)*log(1/epsilon),
        return ((size_t)std::round((double)n/(double)k) * std::log10(1/epsilon)) + 1;
    }

public:
    MaxKCoverEngine(int k) 
    {
        this->k = k;
    };

    ~MaxKCoverEngine() {
        delete this->finder;
    }

    MaxKCoverEngine* useStochasticGreedy(double e)
    {
        this->epsilon = e;
        this->usingStochastic = true;
        return this;
    }

    MaxKCoverEngine* useLazyGreedy(std::unordered_map<int, std::unordered_set<int>>& data)
    {
        this->finder = new LazyGreedy(data);

        return this;
    }

    MaxKCoverEngine* useNaiveGreedy(std::unordered_map<int, std::unordered_set<int>>& data)
    {
        this->finder = new NaiveGreedy(data);

        return this;
    }

    std::pair<std::vector<unsigned int>, int> run_max_k_cover(std::unordered_map<int, std::unordered_set<int>>& data, ssize_t theta)
    {
        std::vector<unsigned int> res(this->k, -1);
        ripples::Bitmask<int> covered(theta);

        size_t subset_size = (this->usingStochastic) ? this->getSubsetSize(data.size(), this->k, this->epsilon) : data.size() ;

        std::vector<unsigned int>* all_vertices = new std::vector<unsigned int>();  
        for (const auto & l : data) { all_vertices->push_back(l.first); }
        this->finder->setSubset(all_vertices, subset_size);

        int uniqueCounted = 0;
        for (int currentSeed = 0; currentSeed < k; currentSeed++)
        {
            if (this->usingStochastic)
            {
                reorganizeVertexSet(all_vertices, subset_size, res);
                this->finder->reloadSubset();
            }

            uniqueCounted += finder->findNextInfluential(
                res, currentSeed, covered, theta
            );
        }
        
        delete all_vertices;
        return std::make_pair(res, uniqueCounted);
    }
};