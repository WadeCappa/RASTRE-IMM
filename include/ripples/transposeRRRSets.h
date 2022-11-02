#include <vector>
#include <unordered_set>
#include <mutex>
#include <iostream>

typedef struct linearizedSetsSize {
    int count;
    std::vector<int> countPerProcess;
} LinearizedSetsSize;

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

    LinearizedSetsSize count(std::vector<int> vertexToProcessor, int p) {
        std::vector<int> countPerProcess(p); 
        LinearizedSetsSize* setSize = new LinearizedSetsSize();
        int count = 0;
        for (int i = 0; i < vertexToProcessor.size(); i++) {
            count += (*sets)[i].second->size() + 2;
            countPerProcess[vertexToProcessor[i]] += (*sets)[i].second->size() + 2;
        }   
        setSize->count = count;
        setSize->countPerProcess = countPerProcess;

        return *setSize;
    }

    std::vector<int> buildPartialSum(std::vector<int> processorSizes) {
        std::vector<int>* partialSum = new std::vector<int>(processorSizes.size(), 0);

        for (int i = 1; i < processorSizes.size(); i++) {
            (*partialSum)[i] = (*partialSum)[i - 1] + processorSizes[i - 1];
        }
        
        return *partialSum;
    }

    int* linearize(std::vector<int> vertexToProcessor, std::vector<int> dataStartPartialSum, int totalData, int p) 
    {
        int* linearTRRRSets = new int[totalData];

        #pragma omp for
        for (int rank = 0; rank < p; rank++) {
            linearize(linearTRRRSets, vertexToProcessor, dataStartPartialSum, totalData, rank);
        }
        return linearTRRRSets;
    }

    void linearize(int* linearTRRRSets, std::vector<int> vertexToProcessor, std::vector<int> dataStartPartialSum, int totalData, int rank ) 
    {
        // TODO: Calculate the partial sum, set the index to the correct starting position for the incoding data. 
        int index = dataStartPartialSum[rank];
        for (int i = 0; i < sets->size(); i++) {
            if (vertexToProcessor[i] == rank) {
                linearTRRRSets[index++] = i;
                for (const auto& RRRid: *((*sets)[i].second)) {
                    linearTRRRSets[index++] = RRRid;
                }
                linearTRRRSets[index++] = -1;
            }
        }
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