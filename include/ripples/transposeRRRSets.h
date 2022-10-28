#include <vector>
#include <unordered_set>
#include <mutex>

template <typename GraphTy>
class TransposeRRRSets
{
    private:
    std::vector<std::pair<std::mutex*, std::unordered_set<typename GraphTy::vertex_type>*>> *sets;

    public:
    TransposeRRRSets(int numberOfVertices) 
    {
        sets = new std::vector<std::pair<std::mutex*, std::unordered_set<typename GraphTy::vertex_type>*>>(numberOfVertices);
        for (int i = 0; i < numberOfVertices; i++) {
            (*sets)[i].first = new std::mutex();
            (*sets)[i].second = new std::unordered_set<typename GraphTy::vertex_type>();
        }
    }

    ~TransposeRRRSets() 
    {
        delete sets;
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