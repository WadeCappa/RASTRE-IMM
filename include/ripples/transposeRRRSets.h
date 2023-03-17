#include <vector>
#include <set>
#include <mutex>
#include <iostream>

template <typename GraphTy>
class TransposeRRRSets
{
    public:
    std::vector<std::pair<std::mutex*, std::set<typename GraphTy::vertex_type>*>> *sets;

    public:
    TransposeRRRSets(int numberOfVertices) 
    {
        sets = new std::vector<std::pair<std::mutex*, std::set<typename GraphTy::vertex_type>*>>(numberOfVertices);
        for (int i = 0; i < numberOfVertices; i++) {
            (*sets)[i].first = new std::mutex();
            (*sets)[i].second = new std::set<typename GraphTy::vertex_type>();
        }
    }

    ~TransposeRRRSets() 
    {
        for (int i = 0; i < sets->size(); i++) {
            free(sets->at(i).first);
            free(sets->at(i).second);
        }

        free(sets);
    }

    void addToSet(int index, typename GraphTy::vertex_type vertex) 
    {
        (*sets)[index].first->lock(); 
        (*sets)[index].second->insert(vertex);
        (*sets)[index].first->unlock();
    }

    auto getBegin() 
    {
        return sets->begin();
    }

    auto getEnd() 
    {
        return sets->end();
    }
};

