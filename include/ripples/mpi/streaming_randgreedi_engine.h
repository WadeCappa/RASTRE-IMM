#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <iostream> 
#include <utility>
#include <chrono>
#include <math.h>
#include "ripples/bitmask.h"


class ThresholdBucket 
{
    private:
    ripples::Bitmask<int> localCovered;
    std::vector<std::pair<int,int>> seeds;
    double marginalGainThreshold;
    int foundSeeds;

    public:
    ThresholdBucket(size_t theta, int deltaZero, int k, double epsilon, size_t numBucket) 
        : localCovered(theta), seeds(k)
    {
        this->marginalGainThreshold = ( deltaZero / ( 2 * k )) * pow(1 + epsilon, numBucket);
    }

    size_t getUtility() 
    {
        return localCovered.popcount();
    }

    std::vector<std::pair<int,int>> getSeeds()
    {
        return this->seeds;
    }

    int getTotalCovered()
    {
        return foundSeeds;
    }

    bool attemptInsert(const std::pair<int, std::unordered_set<int>*>& element) 
    {
        std::unordered_set<int> temp;
        for (const int & e : *(element.second)) {
            if (!localCovered.get(e))
                temp.insert(e);
        }
        if (temp.size() >= marginalGainThreshold) {
            for (const int & e : temp) {
                localCovered.set(e);
            }

            this->seeds[foundSeeds].first = element.first;
            this->seeds[foundSeeds].second = temp.size();
            foundSeeds++;                 
            return true;
        }
        return false;
    }
};


class StreamingRandGreedIEngine 
{
    private:
    std::vector<ThresholdBucket> buckets;
    size_t theta;
    int k;
    int world_size;
    std::vector<std::pair<MPI_Request, int*>> receive_buffers;

    public:
    StreamingRandGreedIEngine(size_t theta, int deltaZero, int k, double epsilon, size_t totalBuckets, int world_size) 
        : receive_buffers(world_size)
    {

        for (size_t i = 0; i < totalBuckets; i++) 
        {
            ThresholdBucket b(theta, deltaZero, k, epsilon, i);
            buckets.push_back(b);
        }

        this->theta = theta;
        this->k = k;
        this->world_size = world_size;

        // initialization step
        for (int i = 0; i < world_size; i++) 
        {
            // In the future this could be a bitmap to save space and lower communciation times
            this->receive_buffers[i].second = new int[this->theta];

            MPI_Irecv (
                this->receive_buffers[i].second,
                this->theta, MPI_INT, i, 0,
                MPI_COMM_WORLD,
                &(this->receive_buffers[i].first)
            );
        }
    }

    ~StreamingRandGreedIEngine() {}

    int getElements(std::vector<std::pair<int, std::unordered_set<int>*>>& result)
    {
        int received_elements = 0;

        for (int i = 0; i < this->world_size; i++) 
        {
            int flag;
            MPI_Status status;
            MPI_Test(&(this->receive_buffers[i].first), &flag, &status);
            if (flag == 1) {
                received_elements++;
                std::unordered_set<int>* received_data = new std::unordered_set<int>();

                for (int e = *(this->receive_buffers[i].second + 1); e != -1; e = *(&e + 1)) 
                {
                    received_data->insert(e);
                }

                result.push_back(std::make_pair(*(this->receive_buffers[i].second), received_data));

                MPI_Irecv (
                    this->receive_buffers[i].second,
                    this->theta, MPI_INT, i, 0, 
                    MPI_COMM_WORLD, 
                    &(this->receive_buffers[i].first)
                );
            }
        }

        return received_elements;
    }

    int insertElements(const std::vector<std::pair<int, std::unordered_set<int>*>>& elements) 
    {
        int totalInsertions = 0;

        #pragma omp parallel for reduction(+:totalInsertions)
        for (size_t t = 0; t < buckets.size(); t++) 
        {
            for (const auto & element : elements) 
            {
                if (buckets[t].getTotalCovered() < this->k && element.first != -1)
                {
                    totalInsertions += int(buckets[t].attemptInsert(element));
                }
            }
        }

        return totalInsertions;
    }

    std::pair<std::vector<unsigned int>, int> stream() 
    {
        // While finished_data_streams < total_streams {
            // 1) Iterate over your vector of receive buffers, std::vector<std::pair<status, int*>> 
            //  for every buffer that has data (test status) read it into another vector of type 
            //  std::vector<std::pair<int, std::unordered_map>> that holds your more permanent data. 

            // 2) While reading stream data, if the tag sent with the data is a 1, this means stop 
            //  reading from this stream. 

            // 3) Using the data read into permanent storage from the last step, in parallel iterate
            //  through all buckets in this.buckets. Loop again over all received data and try to instert
            //  each element into the bucket (which the bucket should independently handle). 

            // 4) Buckets only insert IFF they contain less than k elements and the inserted element
            //  meets this bucket's threshold requirement. 

            // 5) Handle stop conditions. For now the simplist (and least efficent) case will be assumed
            //  that the algorithm only stops when all streams have stopped. 
        // }
 
        // streaming algorithm start
        int received_elements = 0;
        while (received_elements < this->k * world_size) 
        {
            std::vector<std::pair<int, std::unordered_set<int>*>> elements;

            received_elements += this->getElements(elements);

            this->insertElements(elements);

            for (const auto & e : elements)
            {
                delete e.second;
            }
        }

        std::pair<std::vector<unsigned int>, int> result;

        int max_covered = 0;
        int max_covered_index = 0;

        for (int i = 0; i < buckets.size(); i++)
        {
            if (max_covered > this->buckets[i].getUtility()) {
                max_covered = this->buckets[i].getUtility();
                max_covered_index = i;
            }
        }
        
        std::vector<unsigned int>* seeds = new std::vector<unsigned int>();
        for (const auto & p : this->buckets[max_covered_index].getSeeds()) {
            seeds->push_back(p.first);
        }

        return std::make_pair(*seeds, max_covered);
    }

};