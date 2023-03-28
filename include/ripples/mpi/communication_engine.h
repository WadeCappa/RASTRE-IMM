#include <vector>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
#include <mutex>
#include <iostream> 
#include <utility>
#include <chrono>
#include <cmath>

typedef struct linearizedSetsSize {
    int count;
    std::vector<int> countPerProcess;
} LinearizedSetsSize;

template <typename GraphTy>
class CommunicationEngine
{
    private: 

    public:
    CommunicationEngine() {}
    ~CommunicationEngine() {}

    // Returns total count, modifies the countPerProcess vector;
    size_t count(std::vector<unsigned int>& countPerProcess, TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, int p) 
    {
        std::vector<std::mutex> locksPerProcess(p);
        countPerProcess.resize(p, 0);

        size_t count = 0;
        #pragma omp parallel for reduction(+:count)
        for (int i = 0; i < vertexToProcessor.size(); i++) 
        {
            count += tRRRSets.sets[i].second.size() + 2;

            locksPerProcess[vertexToProcessor[i]].lock();
            countPerProcess[vertexToProcessor[i]] += tRRRSets.sets[i].second.size() + 2;
            locksPerProcess[vertexToProcessor[i]].unlock();
        }   

        return count;
    }

    size_t countPerProcessForBatchSend(std::vector<unsigned int>& countPerProcess, TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, int p, size_t block_size) 
    {
        std::vector<std::mutex> locksPerProcess(p);
        countPerProcess.resize(p, 0);

        size_t count = 0;
        #pragma omp parallel for reduction(+:count)
        for (size_t i = 0; i < vertexToProcessor.size(); i++) 
        {
            count += tRRRSets.sets[i].second.size() + 2;

            locksPerProcess[vertexToProcessor[i]].lock();
            countPerProcess[vertexToProcessor[i]] += tRRRSets.sets[i].second.size() + 2;
            locksPerProcess[vertexToProcessor[i]].unlock();
        }   

        size_t block_tax = 0;
        for (auto & processCount : countPerProcess)
        {
            block_tax += (processCount % block_size);
            processCount += (processCount % block_size);
        }

        return count + block_tax;
    }

    void buildPrefixSum(std::vector<unsigned int>& prefixSum, unsigned int* v, int p) {
        prefixSum.resize(p, 0);

        for (int i = 1; i < p; i++) {
            prefixSum[i] = prefixSum[i - 1] + v[i - 1];
        }
    }

    void SendToGlobal(int* data, int size, int tag) {
        MPI_Send(data, size, MPI_INT, 0, tag,
            MPI_COMM_WORLD);
    }

    

    size_t linearizeLocalSeeds(std::vector<unsigned int>& linearAggregateSets, const std::map<int, std::vector<int>> &aggregateSets, const std::vector<unsigned int>& localSeeds, const size_t total_utility) 
    {
        std::vector<std::pair<int, int>> setsPrefixSum;
        size_t runningSum = 0;

        for (const auto & seed : localSeeds) {
            auto seed_set = aggregateSets.at(seed);
            setsPrefixSum.push_back(std::make_pair(seed, runningSum));
            runningSum += seed_set.size() + 2;
        }

        size_t totalData = runningSum + 1;
        linearAggregateSets.resize(totalData);

        #pragma omp parallel for
        for (int setIndex = 0; setIndex < setsPrefixSum.size(); setIndex++) {
            linearAggregateSets[setsPrefixSum[setIndex].second] = setsPrefixSum[setIndex].first;
            int offset = setsPrefixSum[setIndex].second + 1;
            for (const auto & RRRSetID : aggregateSets.at(setsPrefixSum[setIndex].first)) {
                linearAggregateSets[offset++] = RRRSetID;
            }
            linearAggregateSets[setsPrefixSum[setIndex].second + aggregateSets.at(setsPrefixSum[setIndex].first).size() + 1] = -1;
        }

        // mark end of process seeds
        linearAggregateSets[totalData - 2] = -2;

        // mark total utility of local process
        linearAggregateSets[totalData - 1] = total_utility;

        return totalData;
    }

    std::pair<int, int*> LinearizeSingleSeed(const int vertexID, const std::vector<int>& set) 
    {
        int* linearAggregateSets = new int[set.size()+2];
        int offset = 1;

        *linearAggregateSets = vertexID;

        for (const auto & RRRSetID : set) {
            *(linearAggregateSets + offset++) = RRRSetID;
        }

        *(linearAggregateSets + offset) = -1;

        return std::make_pair(set.size()+2, linearAggregateSets);
    }

    void linearize(unsigned int* linearTRRRSets, TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, std::vector<unsigned int> dataStartPartialSum, size_t totalData, int p, size_t block_size) 
    {
        #pragma omp parallel for
        for (int rank = 0; rank < p; rank++) {
            linearizeRank(tRRRSets, linearTRRRSets, vertexToProcessor, dataStartPartialSum, rank, totalData);
        }
    }

    void linearizeRank(TransposeRRRSets<GraphTy> &tRRRSets, unsigned int* linearTRRRSets, std::vector<int> vertexToProcessor, std::vector<unsigned int> dataStartPartialSum, int rank, size_t totalData) 
    {
        size_t index = dataStartPartialSum[rank];
        for (size_t i = 0; i < tRRRSets.sets.size(); i++) {
            if (vertexToProcessor[i] == rank) {
                linearTRRRSets[index++] = i;
                for (const auto& RRRid: tRRRSets.sets[i].second) {
                    linearTRRRSets[index++] = RRRid;
                }
                linearTRRRSets[index++] = -1;
            }
        }

        size_t stop = rank < (dataStartPartialSum.size() - 1) ? dataStartPartialSum[rank+1] : totalData ;
        for (; index < stop; index++) {
            linearTRRRSets[index] = -1;
        }
    }

    void DEBUG_printLinearizedSets(int* linearizedData, int totalData)
    {
        for (int i = 0; i < totalData; i++) {
            std::cout << linearizedData[i] << " ";
            if (linearizedData[i] == -1) {
              std::cout << std::endl;
            }
          }
          std::cout << std::endl;
    }

    /// @brief every process calls this function after they are sent the bulk data from MPI_alltoallV. Each
    ///     process only gets a certain set of values. 
    /// @param tRRRSets is the local process's tRRRSets. This should be a new object (i.e, empty).
    /// @param data is the data collected from MPI_alltoallV
    /// @param receivedDataSizes is the data collected from MPI_alltoall
    /// @param p number of processes
    /// @param RRRIDsPerProcess the upper bound of the maximum number of RRRIDs that each process is responsible for generating
    void aggregateTRRRSets(std::map<int, std::vector<int>> &aggregateSets, unsigned int* data, unsigned int* receivedDataSizes, int p, ssize_t RRRIDsPerProcess)
    {
        size_t totalData = 0;
        for (int i = 0; i < p; i++) {
            totalData += *(receivedDataSizes + i);
        }            

        int receivedDataRank = 0;

        // cycle over data
        int vertexID = *data;   
        aggregateSets.insert({ vertexID, std::vector<int>() });
        for (size_t rankDataProcessed = 1, i = 1; i < totalData - 1; i++, rankDataProcessed++) {
            if (*(data + i) == -1 && *(data + i-1) != -1) {
                vertexID = *(data + ++i);
                rankDataProcessed++;
                aggregateSets.insert({ vertexID, std::vector<int>() });
            }

            else if (*(data + i) != -1) {
                // add each RRRSetID to the map indexed by the target vertex id. 
                // The RRRSetIDs need to be modified such that they are unique from the sent process. 
                // This is possible because the order is known and the expected sizes of data that should be sent. 
                aggregateSets[vertexID].push_back(*(data + i) + (RRRIDsPerProcess * receivedDataRank));
            }

            // track the current receiving process rank
            if (rankDataProcessed > *(receivedDataSizes + receivedDataRank)) {
                receivedDataRank += 1;
                rankDataProcessed = 1;
            }
        }
    }

    std::vector<std::pair<unsigned int, std::vector<unsigned int>*>>* aggregateLocalKSeeds(std::map<int, std::vector<int>> &bestMKSeeds, unsigned int* data, size_t totalData)
    {
        // tracks total utility of each local process
        std::vector<std::pair<unsigned int, std::vector<unsigned int>*>>* local_seeds = new std::vector<std::pair<unsigned int, std::vector<unsigned int>*>>();
        std::vector<unsigned int>* current_seeds = new std::vector<unsigned int>();

        // cycle over data
        int vertexID = *data;   
        current_seeds->push_back(vertexID);
        bestMKSeeds.insert({ vertexID, std::vector<int>() });

        for (size_t rankDataProcessed = 1, i = 1; i < totalData - 1; i++) {
            if (*(data + i) == -2) {
                local_seeds->push_back(std::make_pair(*(data+i+1), current_seeds));
                current_seeds = new std::vector<unsigned int>();
                i++;
                vertexID = *(data + ++i);
                current_seeds->push_back(vertexID);
                bestMKSeeds.insert({ vertexID, std::vector<int>() });
            }

            else if (*(data + i) == -1) {
                current_seeds->push_back(vertexID);
                vertexID = *(data + ++i);
                bestMKSeeds.insert({ vertexID, std::vector<int>() });
            }

            else {
                bestMKSeeds[vertexID].push_back(*(data + i));
            }
        }

        return local_seeds;
    }


    void getProcessSpecificVertexRRRSets(
        std::map<int, std::vector<int>> &aggregateSets, 
        unsigned int* linearizedData,
        unsigned int* countPerProcess,
        int p,
        ssize_t RRRIDsPerProcess,
        MPI_Datatype& batch_int,
        int block_size
    ) {
        unsigned int* receiveSizes = new unsigned int[p];

        MPI_Alltoall(countPerProcess, 1, MPI_UNSIGNED, receiveSizes, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

        size_t totalReceivingBlocks = 0;
        for (int i = 0; i < p; i++) {
            totalReceivingBlocks += (unsigned int)*(receiveSizes + i);
        }
        
        unsigned int* linearizedLocalData = new unsigned int[totalReceivingBlocks * block_size];

        std::vector<unsigned int> sendPrefixSum;
        buildPrefixSum(sendPrefixSum, countPerProcess, p);
        std::vector<unsigned int> receivePrefixSum; 
        buildPrefixSum(receivePrefixSum, receiveSizes, p);

        // call alltoall_v
        MPI_Alltoallv(
            linearizedData, 
            (int*)countPerProcess, 
            (int*)sendPrefixSum.data(), 
            batch_int, 
            linearizedLocalData, 
            (int*)receiveSizes, 
            (int*)receivePrefixSum.data(),
            batch_int,
            MPI_COMM_WORLD
        );

        for (int i = 0; i < p; i++) {
            receiveSizes[i] *= block_size;
        }

        aggregateTRRRSets(aggregateSets, linearizedLocalData, receiveSizes, p, RRRIDsPerProcess);
        
        delete receiveSizes;
        delete linearizedLocalData;
    }

    size_t aggregateAggregateSets(std::vector<unsigned int>& aggregatedSeeds, size_t totalGatherData, int world_size, unsigned int* localLinearAggregateSets)
    {
        unsigned int* gatherSizes = new unsigned int[world_size];
        MPI_Allgather(
            &totalGatherData, 1, MPI_INT,
            gatherSizes, 1, MPI_INT, MPI_COMM_WORLD
        );
        
        // mpi_gatherv for all buffers of linearized aggregateSets. 
        size_t totalData = 0;
        for (int i = 0; i < world_size; i++) {
            totalData += (unsigned int)gatherSizes[i];
        }  

        std::vector<unsigned int> displacements;
        buildPrefixSum(displacements, gatherSizes, world_size);

        aggregatedSeeds.resize(totalData);

        MPI_Gatherv(
            localLinearAggregateSets,
            totalGatherData,
            MPI_UNSIGNED,
            (int*)aggregatedSeeds.data(),
            (int*)gatherSizes,
            (int*)displacements.data(),
            MPI_UNSIGNED,
            0,
            MPI_COMM_WORLD
        );

        delete gatherSizes;
        return totalData;
    }

    void distributeF(double* f)
    {
        MPI_Bcast(f, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
};


