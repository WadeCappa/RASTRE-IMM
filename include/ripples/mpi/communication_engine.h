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
    const int block_size = (1024 / sizeof(unsigned int)); // 256 ints per block

    const unsigned int world_size;
    const unsigned int world_rank;

    public:
    CommunicationEngine
    (
        const unsigned int world_size, 
        const unsigned int world_rank
    ) : world_size(world_size), world_rank(world_rank)
    {
        // TODO: Extract all MPI_Datatype operations into constructor and deconstructor. Last time you tried
        //  this you got MPI_DATATYPE_NULL errors from mpich. Figure out how to prevent this.
    }

    ~CommunicationEngine() 
    {
    }

    size_t GetSendReceiveBufferSize(const size_t data_size) const
    {
        return (size_t)(data_size + (size_t)(data_size % this->block_size == 0 ? 0 : (this->block_size - (data_size % this->block_size))));
    }

    unsigned int GetRank() const
    {
        return this->world_rank;
    }

    unsigned int GetSize() const
    {
        return this->world_size;
    }

    int TestSend(MPI_Request &request) const
    {
        MPI_Status status;
        int flag = 0;
        MPI_Test(&request, &flag, &status);
        return flag;
    }

    void SendNextSeed(
        const unsigned int* data, 
        const size_t data_size, 
        const unsigned int tag
    ) const
    {
        MPI_Datatype batch_int;
        MPI_Type_contiguous( this->block_size, MPI_INT, &(batch_int) );
        MPI_Type_commit(&(batch_int));

        MPI_Send (
            data,
            data_size / this->block_size,
            batch_int, 0,
            tag,
            MPI_COMM_WORLD
        );
        MPI_Type_free(&(batch_int));
    }

    void QueueNextSeed(
        const unsigned int* data, 
        const size_t data_size, 
        const unsigned int tag, 
        MPI_Request &request
    ) const
    {
        MPI_Datatype batch_int;
        MPI_Type_contiguous( this->block_size, MPI_INT, &(batch_int) );
        MPI_Type_commit(&(batch_int));
        MPI_Isend (
            data,
            data_size / this->block_size,
            batch_int, 
            0,
            tag,
            MPI_COMM_WORLD,
            &request
        );
        MPI_Type_free(&(batch_int));
    }

    size_t GetBlockSize() const
    {
        return this->block_size;
    }

    std::vector<int> DistributeVertices (
        const bool streaming,
        const GraphTy &G
    ) const
    {
        size_t nodes = G.num_nodes();

        std::vector<int> vertex_mapping(nodes, -1);
        unsigned int seed = (unsigned int)time(0);
        MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        std::uniform_int_distribution<int> uniform_distribution(streaming ? 1 : 0, this->world_size - 1);
        std::default_random_engine number_selecter(seed);

        for (int i = 0; i < nodes; i++) {
            vertex_mapping[i] = uniform_distribution(number_selecter);
        }

        return vertex_mapping;
    }

    void AggregateThroughAllToAll(
        const TransposeRRRSets<GraphTy> &tRRRSets,
        const std::vector<int> &vertexToProcess,
        std::vector<size_t> &old_sampling_sizes,
        const size_t localThetaPrime,
        std::map<int, std::vector<int>> &aggregateSets
    ) const
    {
        std::vector<unsigned int> countPerProcess;
        size_t totalCount = this->countPerProcessForBatchSend(countPerProcess, tRRRSets, vertexToProcess, old_sampling_sizes);

        std::vector<unsigned int> psum;
        this->buildPrefixSum(psum, countPerProcess.data(), this->world_size);

        std::vector<unsigned int> linear_sets = this->LinearizeColumnsInRRRSetMatrix(
            tRRRSets, 
            vertexToProcess, 
            psum, 
            totalCount, 
            old_sampling_sizes
        );

        for (size_t i = 0; i < old_sampling_sizes.size(); i++)
        {
            old_sampling_sizes[i] = tRRRSets.sets[i].second.size();
        }

        // std::cout << "linearizing data, rank = " << world_rank << std::endl;

        for (auto & processCount : countPerProcess) {
            processCount = processCount / block_size;
        }
        
        this->GetProcessSpecificVertexRRRSets(aggregateSets, linear_sets.data(), countPerProcess.data(), localThetaPrime);
    }

    // void BuildBufferForStreamingSend(std::vector<unsigned int>& buffer) const
    // {   
    //     size_t total_data = this->GetSendReceiveBufferSize((size_t)sendData.first);

    //     buffer.resize(total_data);

    //     size_t i = 0;

    //     for (; i < sendData.first; i++)
    //     {
    //         buffer[i] = sendData.second[i];
    //     }

    //     for (; i < total_data; i++)
    //     {
    //         buffer[i] = -1; //TODO: Verify that this end of message flag is correct
    //     }

    //     delete sendData.second;
    // }

    void QueueReceive(unsigned int* buffer, size_t max_receive_size, MPI_Request &request) const 
    {
        MPI_Datatype batch_int;
        MPI_Type_contiguous( this->block_size, MPI_INT, &(batch_int) );
        MPI_Type_commit(&(batch_int));
        MPI_Irecv(
            buffer,
            max_receive_size / this->block_size,
            batch_int,
            MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            MPI_COMM_WORLD,
            &request
        );
        MPI_Type_free(&(batch_int));
    }

    void WaitForSend(MPI_Request &request) const 
    {
        MPI_Status status;
        MPI_Wait(&request, &status);
    }

    size_t GatherPartialSolutions(
        const std::pair<std::vector<unsigned int>, ssize_t> &localSeeds,
        const std::map<int, std::vector<int>> &aggregateSets,
        std::vector<unsigned int> &globalAggregation,
        const MPI_Comm &commWorld
    ) const
    {
        std::vector<unsigned int> linearAggregateSets;
        size_t totalLinearLocalSeedsData = this->linearizeLocalSeeds(linearAggregateSets, aggregateSets, localSeeds.first, localSeeds.second);
        return this->aggregateAggregateSets(globalAggregation, totalLinearLocalSeedsData, linearAggregateSets.data(), commWorld);
    }

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
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        const std::vector<int> vertexToProcessor, 
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
            size_t blocks_to_add = this->block_size - (processCount % this->block_size);
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
        const size_t total_utility
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
        size_t total_block_data = totalData + (this->block_size - (totalData % this->block_size));

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

    std::pair<size_t, int*> LinearizeSingleSeed(const int vertexID, const std::vector<int>& set) const
    {
        int* linearAggregateSets = new int[set.size()+2];
        size_t offset = 1;

        *linearAggregateSets = vertexID;

        for (const auto & RRRSetID : set) {
            *(linearAggregateSets + offset++) = RRRSetID;
        }

        *(linearAggregateSets + offset) = -1;

        return std::make_pair(set.size()+2, linearAggregateSets);
    }

    std::vector<unsigned int> LinearizeColumnsInRRRSetMatrix(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        const std::vector<int> &vertexToProcessor, 
        const std::vector<unsigned int> &dataStartPartialSum, 
        const size_t total_data, 
        const std::vector<size_t>& old_sizes
    ) const
    {
        std::vector<unsigned int> linear_sets(total_data);

        #pragma omp parallel for
        for (int rank = 0; rank < this->world_size; rank++) {
            this->LinearizeRank(tRRRSets, linear_sets, vertexToProcessor, dataStartPartialSum, rank, old_sizes);
        }

        return linear_sets;
    }

    void LinearizeRank(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        std::vector<unsigned int>& linearTRRRSets, 
        const std::vector<int> &vertexToProcessor, 
        const std::vector<unsigned int> &dataStartPartialSum, 
        const int rank, 
        const std::vector<size_t>& old_sizes
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

        size_t stop = rank < (dataStartPartialSum.size() - 1) ? dataStartPartialSum[rank+1] : linearTRRRSets.size() ;
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
    void AggregateTRRRSets(
        std::map<int, std::vector<int>> &aggregateSets, 
        unsigned int* data, 
        unsigned int* receivedDataSizes, 
        size_t RRRIDsPerProcess
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

    std::vector<std::pair<unsigned int, std::vector<unsigned int>>> aggregateLocalKSeeds (
        std::map<int, std::vector<int>> &bestMKSeeds, 
        unsigned int* data, 
        size_t totalData
    ) const
    {
        // tracks total utility of each local process
        std::vector<std::pair<unsigned int, std::vector<unsigned int>>> local_seeds;
        std::vector<unsigned int> current_seeds;

        // cycle over data
        int vertexID = *data;   
        current_seeds.push_back(vertexID);
        bestMKSeeds.insert({ vertexID, std::vector<int>() });

        for (size_t i = 1; i < totalData - 1; i++) {

            // std::cout << "processing byte "<< i << " out of " << totalData << std::endl;
            if (*(data + i) == -2) {
                local_seeds.push_back(std::make_pair(*(data + ++i), std::vector<unsigned int>(current_seeds)));
                current_seeds.empty();

                i++;
                while (*(data + i) == -3 && i < totalData)
                {
                    i++;
                    if (i >= totalData) {
                        return local_seeds;
                    }
                }
                
                vertexID = *(data + i);
                current_seeds.push_back(vertexID);
                bestMKSeeds.insert({ vertexID, std::vector<int>() });
            }

            else if (*(data + i) == -1) {
                current_seeds.push_back(vertexID);
                vertexID = *(data + ++i);
                bestMKSeeds.insert({ vertexID, std::vector<int>() });
            }

            else {
                bestMKSeeds[vertexID].push_back(*(data + i));
            }
        }

        return local_seeds;
    }


    void GetProcessSpecificVertexRRRSets(
        std::map<int, std::vector<int>> &aggregateSets, 
        unsigned int* linearizedData,
        unsigned int* countPerProcess,
        size_t RRRIDsPerProcess
    ) const {
        unsigned int* receiveSizes = new unsigned int[this->world_size];

        MPI_Alltoall(countPerProcess, 1, MPI_UNSIGNED, receiveSizes, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

        size_t totalReceivingBlocks = 0;
        for (int i = 0; i < this->world_size; i++) {
            totalReceivingBlocks += (unsigned int)*(receiveSizes + i);
	    // std::cout << "rank " << GetRank() << " receiving " << (unsigned int)*(receiveSizes + i) << " blocks from " << i << std::endl;
        }
        
        unsigned int* linearizedLocalData = new unsigned int[totalReceivingBlocks * this->block_size];

        std::vector<unsigned int> sendPrefixSum;
        buildPrefixSum(sendPrefixSum, countPerProcess, this->world_size);
        std::vector<unsigned int> receivePrefixSum; 
        buildPrefixSum(receivePrefixSum, receiveSizes, this->world_size);

        // call alltoall_v
        MPI_Datatype batch_int;
        MPI_Type_contiguous( this->block_size, MPI_INT, &(batch_int) );
        MPI_Type_commit(&(batch_int));
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
        MPI_Type_free(&(batch_int));

        for (int i = 0; i < this->world_size; i++) {
            receiveSizes[i] *= this->block_size;
        }

        // TODO: make this not terrible
	//  for streaming this might crash so you may need to add a dummy value to aggregate sets.
        //if (this->world_rank == 0)
        //{
        //    aggregateSets.insert({0, std::vector<int>()});
        //}   

    	this->AggregateTRRRSets(aggregateSets, linearizedLocalData, receiveSizes, RRRIDsPerProcess);
        
        delete receiveSizes;
        delete linearizedLocalData;
    }

    size_t aggregateAggregateSets(
        std::vector<unsigned int>& aggregatedSeeds, 
        size_t totalGatherData, 
        unsigned int* localLinearAggregateSets,
        const MPI_Comm &commWorld
    ) const
    {
        unsigned int* gatherSizes = new unsigned int[this->world_size];
        gatherSizes[0] = -1;

        size_t total_block_send = totalGatherData / this->block_size;
        // std::cout << "totalGatherData " << totalGatherData << std::endl; 

        MPI_Gather(
            &total_block_send, 1, MPI_INT,
            gatherSizes, 1, MPI_INT, 0, commWorld
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

        aggregatedSeeds.resize(totalData * this->block_size);

        // for the leveled implementation of this communication, we may have to use alltoallv 
        //  and formulate the communication as a matrix manipulation. Look more into this, create
        //  some basic tests scripts.
        MPI_Datatype batch_int;
        MPI_Type_contiguous( this->block_size, MPI_INT, &(batch_int) );
        MPI_Type_commit(&(batch_int));
        MPI_Gatherv(
            localLinearAggregateSets,
            total_block_send,
            batch_int,
            (int*)aggregatedSeeds.data(),
            (int*)gatherSizes,
            (int*)displacements.data(),
            batch_int,
            0,
            commWorld
        );
        MPI_Type_free(&(batch_int));

        delete gatherSizes;
        return totalData * this->block_size;
    }

    void distributeF(double* f) const
    {
        MPI_Bcast(f, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
};

template <typename GraphTy>
class CommunicationEngineBuilder
{
    public:
    static CommunicationEngine<GraphTy> BuildCommunicationEngine()
    {
        int world_size, world_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

        return CommunicationEngine<GraphTy>(world_size, world_rank);
    }
};
