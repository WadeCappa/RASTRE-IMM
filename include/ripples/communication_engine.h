#include <vector>
#include <unordered_set>
#include <mutex>
#include <iostream>

typedef struct linearizedSetsSize {
    int count;
    std::vector<int> countPerProcess;
} LinearizedSetsSize;

template <typename GraphTy>
class CommunicationEngine
{
    public:
    CommunicationEngine() {}
    ~CommunicationEngine() {}

    LinearizedSetsSize count(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, int p) {
        std::vector<int> countPerProcess = *(new std::vector<int>(p)); 
        LinearizedSetsSize* setSize = new LinearizedSetsSize();
        int count = 0;
        for (int i = 0; i < vertexToProcessor.size(); i++) {
            count += (*tRRRSets.sets)[i].second->size() + 2;
            countPerProcess[vertexToProcessor[i]] += (*tRRRSets.sets)[i].second->size() + 2;
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

    int* linearize(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, std::vector<int> dataStartPartialSum, int totalData, int p) 
    {
        int* linearTRRRSets = new int[totalData];

        #pragma omp for
        for (int rank = 0; rank < p; rank++) {
            linearize(tRRRSets, linearTRRRSets, vertexToProcessor, dataStartPartialSum, totalData, rank);
        }
        return linearTRRRSets;
    }

    void linearize(TransposeRRRSets<GraphTy> &tRRRSets, int* linearTRRRSets, std::vector<int> vertexToProcessor, std::vector<int> dataStartPartialSum, int totalData, int rank ) 
    {
        // TODO: Calculate the partial sum, set the index to the correct starting position for the incoding data. 
        int index = dataStartPartialSum[rank];
        for (int i = 0; i < tRRRSets.sets->size(); i++) {
            if (vertexToProcessor[i] == rank) {
                linearTRRRSets[index++] = i;
                for (const auto& RRRid: *((*tRRRSets.sets)[i].second)) {
                    linearTRRRSets[index++] = RRRid;
                }
                linearTRRRSets[index++] = -1;
            }
        }
    }


    void DEBUG_printLinearizedSets(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcesses, int p)
    {
        LinearizedSetsSize setSize = count(tRRRSets, vertexToProcesses, p);
        int* linearizedData = linearize(tRRRSets, vertexToProcesses, buildPartialSum(setSize.countPerProcess), setSize.count, p);
        for (int i = 0; i < setSize.count; i++) {
            std::cout << linearizedData[i] << " ";
            if (linearizedData[i] == -1) {
              std::cout << std::endl;
            }
          }
          std::cout << std::endl;
    }


    int* linearizeTRRRSets(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcesses, int p) 
    {
        LinearizedSetsSize setSize = count(tRRRSets, vertexToProcesses, p);
        return linearize(tRRRSets, vertexToProcesses, buildPartialSum(setSize.countPerProcess), setSize.count, p);
    }
};