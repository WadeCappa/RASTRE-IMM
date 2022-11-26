#include <vector>
#include <unordered_set>
#include <unordered_map>
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

    LinearizedSetsSize* count(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcessor, int p) {
        std::vector<int> countPerProcess = *(new std::vector<int>(p)); 
        LinearizedSetsSize* setSize = new LinearizedSetsSize();
        int count = 0;
        for (int i = 0; i < vertexToProcessor.size(); i++) {
            count += (*tRRRSets.sets)[i].second->size() + 2;
            countPerProcess[vertexToProcessor[i]] += (*tRRRSets.sets)[i].second->size() + 2;
        }   
        setSize->count = count;
        setSize->countPerProcess = countPerProcess;

        return setSize;
    }

    std::vector<int>* buildPrefixSum(int* v, int p) {
        std::vector<int>* prefixSum = new std::vector<int>(p, 0);

        for (int i = 1; i < p; i++) {
            (*prefixSum)[i] = (*prefixSum)[i - 1] + v[i - 1];
        }
        
        return prefixSum;
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

    int* linearizeTRRRSets(TransposeRRRSets<GraphTy> &tRRRSets, std::vector<int> vertexToProcesses, int p) 
    {
        LinearizedSetsSize* setSize = count(tRRRSets, vertexToProcesses, p);
        return linearize(tRRRSets, vertexToProcesses, buildPrefixSum(setSize->countPerProcess), setSize->count, p);
    }

    /// @brief every process calls this function after they are sent the bulk data from MPI_alltoallV. Each
    ///     process only gets a certain set of values. 
    /// @param tRRRSets is the local process's tRRRSets. This should be a new object (i.e, empty).
    /// @param data is the data collected from MPI_alltoallV
    /// @param receivedDataSizes is the data collected from MPI_alltoall
    /// @param p number of processes
    /// @param RRRIDsPerProcess the upper bound of the maximum number of RRRIDs that each process is responsible for generating
    void aggregateTRRRSets(std::unordered_map<int, std::unordered_set<int>> &aggregateSets, int* data, int* receivedDataSizes, int p, int RRRIDsPerProcess)
    {
        int totalData = 0;
        for (int i = 0; i < p; i++) {
            totalData += *(receivedDataSizes + i);
        }            

        int receivedDataRank = 0;

        // cycle over data
        int vertexID = *data;   
        aggregateSets.insert({ vertexID, *(new std::unordered_set<int>()) });
        for (int rankDataProcessed = 1, i = 1; i < totalData - 1; i++, rankDataProcessed++) {
            if (*(data + i) == -1) {
                vertexID = *(data + ++i);
                rankDataProcessed++;
                aggregateSets.insert({ vertexID, *(new std::unordered_set<int>()) });
            }

            else {
                // add each RRRSetID to the map indexed by the target vertex id. 
                // The RRRSetIDs need to be modified such that they are unique from the sent process. 
                // This is possible because the order is known and the expected sizes of data that should be sent. 
                aggregateSets[vertexID].insert(*(data + i) + (RRRIDsPerProcess * receivedDataRank));
            }

            // track the current receiving process rank
            if (rankDataProcessed > *(receivedDataSizes + receivedDataRank)) {
                receivedDataRank += 1;
                rankDataProcessed = 1;
            }
        }
    }


    // All To All
    // 1. World size send buffer containing the value to be sent to every other process
    // 2. World size receive buffer 

    // Send
    // 1. What am I sending
    // 2. How many am I sending
    // 3. Where to count from
    // Receive
    // 4. Where am I receiving
    // 5. How many am I receiving
    // 6. Where to start counting from
    int* allToAll(int* linearizedData, std::vector<int> countPerProcess, int p)
    {
        int* receiveSizes = new int[p];
        MPI_Alltoall(&countPerProcess[0], 1, MPI_INT, receiveSizes, 1, MPI_INT, MPI_COMM_WORLD);

        int totalReceiveSize = 0;
        for (int i = 0; i < p; i++) {
            totalReceiveSize += *(receiveSizes + i);
        }
        int* receiveBuffer = new int[totalReceiveSize];

        // call alltoall_v
        MPI_Alltoallv(
            linearizedData, 
            &countPerProcess[0], 
            buildPrefixSum(countPerProcess), 
            MPI_INT, 
            receiveBuffer, 
            receiveSizes, 
            buildPrefixSum(totalReceiveSize),
            MPI_INT,
            MPI_COMM_WORLD
        );

        return receiveBuffer;
    }

    void getProcessSpecificVertexRRRSets(
        std::unordered_map<int, std::unordered_set<int>> &aggregateSets, 
        int* linearizedData,
        int* countPerProcess,
        int p,
        int RRRIDsPerProcess
    ) {
        int* receiveSizes = new int[p];
        MPI_Alltoall(countPerProcess, 1, MPI_INT, receiveSizes, 1, MPI_INT, MPI_COMM_WORLD);

        int totalReceiveSize = 0;
        for (int i = 0; i < p; i++) {
            totalReceiveSize += *(receiveSizes + i);
        }
        int* linearizedLocalData = new int[totalReceiveSize];

        std::vector<int>* sendPrefixSum = buildPrefixSum(countPerProcess, p);
        std::vector<int>* receivePrefixSum = buildPrefixSum(receiveSizes, p);

        // call alltoall_v
        MPI_Alltoallv(
            linearizedData, 
            &countPerProcess[0], 
            sendPrefixSum->data(), 
            MPI_INT, 
            linearizedLocalData, 
            receiveSizes, 
            receivePrefixSum->data(),
            MPI_INT,
            MPI_COMM_WORLD
        );

        free(sendPrefixSum);
        free(receivePrefixSum);

        aggregateTRRRSets(aggregateSets, linearizedLocalData, receiveSizes, p, RRRIDsPerProcess);
    }

    // TransposeRRRSets<GraphTy> delinearizeToTRRRSets(int* data)
    // {
    //     // build a new tRRRSets object
    //     // cycle over data
    //         // add each RRRSetID to the map indexed by the target vertex id. 
    //         // The RRRSetIDs need to be modified such that they are unique per process. 
    // }
};