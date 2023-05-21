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
    const unsigned int world_size;
    const unsigned int world_rank;

    public:
    CommunicationEngine
    (
        const unsigned int world_size, 
        const unsigned int world_rank
    ) : world_size(world_size), world_rank(world_rank)
    {}

    ~CommunicationEngine() {}

    // Returns total count, modifies the countPerProcess vector;
    size_t count(std::vector<unsigned int>& countPerProcess, TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor) const
    {
        std::vector<std::mutex> locksPerProcess(this->world_size);
        countPerProcess.resize(this->world_size, 0);

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

    size_t countPerProcessForBatchSend(
        std::vector<unsigned int>& countPerProcess, 
        TransposeRRRSets<GraphTy> &tRRRSets, 
        std::vector<int> vertexToProcessor, 
        size_t block_size,
        const std::vector<size_t>& old_sizes
    ) const 
    {
        std::vector<std::mutex> locksPerProcess(this->world_size);
        countPerProcess.resize(this->world_size, 0);

        size_t count = 0;

        #pragma omp parallel for reduction(+:count)
        for (size_t i = 0; i < vertexToProcessor.size(); i++) 
        {
            size_t old_count = old_sizes[i];
            size_t new_size = tRRRSets.sets[i].second.size();
            count += new_size - old_count + 2;
            // count += tRRRSets.sets[i].second.size() + 2;

            locksPerProcess[vertexToProcessor[i]].lock();
            countPerProcess[vertexToProcessor[i]] += new_size - old_count + 2;
            locksPerProcess[vertexToProcessor[i]].unlock();
        }   

        size_t block_tax = 0;
        for (auto & processCount : countPerProcess)
        {
            size_t blocks_to_add = block_size - (processCount % block_size);
            // size_t blocks_to_add = processCount % block_size;
            block_tax += blocks_to_add;
            processCount += blocks_to_add;
        }

        return count + block_tax;
    }

    void buildPrefixSum(std::vector<unsigned int>& prefixSum, unsigned int* v, unsigned int total_count) const
    {
        prefixSum.resize(total_count, 0);

        for (int i = 1; i < total_count; i++) {
            prefixSum[i] = prefixSum[i - 1] + v[i - 1];
        }
    }

    void SendToGlobal(int* data, int size, int tag) const
    {
        MPI_Send(data, size, MPI_INT, 0, tag,
            MPI_COMM_WORLD);
    }

    
    size_t linearizeLocalSeeds(
        std::vector<unsigned int>& linearAggregateSets, 
        const std::map<int, std::vector<int>> &aggregateSets, 
        const std::vector<unsigned int>& localSeeds, 
        const size_t total_utility,
        const int block_size
    ) const
    {
        std::vector<std::pair<int, int>> setsPrefixSum;

        size_t runningSum = 0;
        for (const auto & seed : localSeeds) {
            auto seed_set = aggregateSets.at(seed);
            setsPrefixSum.push_back(std::make_pair(seed, runningSum));
            runningSum += seed_set.size() + 2;
        }

        size_t totalData = runningSum + 1;
        size_t total_block_data = totalData + (block_size - (totalData % block_size));

        linearAggregateSets.resize(total_block_data);

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

        for (int i = totalData; i < total_block_data; i++)
        {
            linearAggregateSets[i] = -3;
        }

        return total_block_data;
    }

    std::pair<int, int*> LinearizeSingleSeed(const int vertexID, const std::vector<int>& set) const
    {
        int* linearAggregateSets = new int[set.size()+2];
        int offset = 1;

        *linearAggregateSets = vertexID;

        for (const auto & RRRSetID : set) {
            *(linearAggregateSets + offset++) = RRRSetID;
        }

        *(linearAggregateSets + offset) = -1;

        std::cout << "sending element " << vertexID << " of size " << set.size()+2 << std::endl;

        // for (int* itr = linearAggregateSets; *itr != -1; itr++)
        // {
        //     std::cout << *itr << ", ";
        // }
        // std::cout << -1 << std::endl;

        return std::make_pair(set.size()+2, linearAggregateSets);
    }

    void linearize(
        unsigned int* linearTRRRSets, 
        TransposeRRRSets<GraphTy> &tRRRSets, 
        std::vector<int> vertexToProcessor, 
        std::vector<unsigned int> dataStartPartialSum, 
        size_t totalData, 
        size_t block_size,
        std::vector<size_t>& old_sizes
    ) const
    {
        #pragma omp parallel for
        for (int rank = 0; rank < this->world_size; rank++) {
            linearizeRank(tRRRSets, linearTRRRSets, vertexToProcessor, dataStartPartialSum, rank, totalData, old_sizes);
        }
    }

    void linearizeRank(
        TransposeRRRSets<GraphTy> &tRRRSets, 
        unsigned int* linearTRRRSets, 
        std::vector<int> vertexToProcessor, 
        std::vector<unsigned int> dataStartPartialSum, 
        int rank, 
        size_t totalData,
        std::vector<size_t>& old_sizes
    ) const
    {
        size_t index = dataStartPartialSum[rank];
        for (size_t i = 0; i < tRRRSets.sets.size(); i++) {
            if (vertexToProcessor[i] == rank) {
                linearTRRRSets[index++] = i;
                for (auto itr = tRRRSets.sets[i].second.begin() + old_sizes[i]; itr != tRRRSets.sets[i].second.end(); itr++) {
                    const auto & RRRid = *itr;
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

    void DEBUG_printLinearizedSets(int* linearizedData, int totalData) const
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
    void aggregateTRRRSets(
        std::map<int, std::vector<int>> &aggregateSets, 
        unsigned int* data, 
        unsigned int* receivedDataSizes, 
        ssize_t RRRIDsPerProcess
    ) const
    {
        size_t totalData = 0;
        for (int i = 0; i < this->world_size; i++) {
            totalData += *(receivedDataSizes + i);
        }            

        int receivedDataRank = 0;

        // cycle over data
        int vertexID = *data; // TODO: Fix all these datatypes, need to be vector_type, will not scale with int as number of nodes increases
        // aggregateSets.insert({ vertexID, std::vector<int>() });
        for (size_t rankDataProcessed = 1, i = 1; i < totalData - 1; i++, rankDataProcessed++) {
            if (*(data + i) == -1 && *(data + i + 1) != -1) {
                vertexID = *(data + ++i);
                rankDataProcessed++;
                // aggregateSets.insert({ vertexID, std::vector<int>() });
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

    std::vector<std::pair<unsigned int, std::vector<unsigned int>*>>* aggregateLocalKSeeds(std::map<int, std::vector<int>> &bestMKSeeds, unsigned int* data, size_t totalData) const
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
                local_seeds->push_back(std::make_pair(*(data + ++i), current_seeds));
                current_seeds = new std::vector<unsigned int>();

                i++;
                while (*(data + i) == -3)
                {
                    i++;
                }
                
                vertexID = *(data + i);
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
        ssize_t RRRIDsPerProcess,
        MPI_Datatype& batch_int,
        int block_size
    ) const {
        unsigned int* receiveSizes = new unsigned int[this->world_size];

        MPI_Alltoall(countPerProcess, 1, MPI_UNSIGNED, receiveSizes, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

        size_t totalReceivingBlocks = 0;
        for (int i = 0; i < this->world_size; i++) {
            totalReceivingBlocks += (unsigned int)*(receiveSizes + i);
        }
        
        unsigned int* linearizedLocalData = new unsigned int[totalReceivingBlocks * block_size];

        std::vector<unsigned int> sendPrefixSum;
        buildPrefixSum(sendPrefixSum, countPerProcess, this->world_size);
        std::vector<unsigned int> receivePrefixSum; 
        buildPrefixSum(receivePrefixSum, receiveSizes, this->world_size);

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

        for (int i = 0; i < this->world_size; i++) {
            receiveSizes[i] *= block_size;
        }

        // TODO: make this not terrible
        if (this->world_rank == 0)
        {
            aggregateSets.insert({0, std::vector<int>()});
        }   
        else 
        {
            aggregateTRRRSets(aggregateSets, linearizedLocalData, receiveSizes, RRRIDsPerProcess);
        }
        
        delete receiveSizes;
        delete linearizedLocalData;
    }

    size_t aggregateAggregateSets(
        std::vector<unsigned int>& aggregatedSeeds, 
        size_t totalGatherData, 
        unsigned int* localLinearAggregateSets,
        int block_size,
        MPI_Datatype& batch_int
    ) const
    {
        unsigned int* gatherSizes = new unsigned int[this->world_size];
        gatherSizes[0] = -1;

        std::cout << "remainder: " << totalGatherData % block_size << std::endl; 
        size_t total_block_send = totalGatherData / block_size;

        MPI_Gather(
            &total_block_send, 1, MPI_INT,
            gatherSizes, 1, MPI_INT, 0, MPI_COMM_WORLD
        );

        std::vector<unsigned int> displacements;
        
        // mpi_gatherv for all buffers of linearized aggregateSets. 
        size_t totalData = 0;
        if (gatherSizes[0] == -1) {
            totalData = 1;
        }
        else{
            for (int i = 0; i < this->world_size; i++) {
                totalData += (unsigned int)gatherSizes[i];
            }  

            buildPrefixSum(displacements, gatherSizes, this->world_size);
        }

        aggregatedSeeds.resize(totalData * block_size);

        MPI_Gatherv(
            localLinearAggregateSets,
            total_block_send,
            batch_int,
            (int*)aggregatedSeeds.data(),
            (int*)gatherSizes,
            (int*)displacements.data(),
            batch_int,
            0,
            MPI_COMM_WORLD
        );

        delete gatherSizes;
        return totalData * block_size;
    }

    void distributeF(double* f) const
    {
        MPI_Bcast(f, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
};


