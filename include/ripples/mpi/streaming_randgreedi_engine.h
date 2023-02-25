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
    std::vector<std::pair<int,int>> seedSetAndMarginalGain;
    double marginalGainThreshold;
    int foundSeeds;

    public:
    ThresholdBucket(size_t theta, int deltaZero, int k, double epsilon, size_t numBucket) 
        : localCovered(theta), seedSetAndMarginalGain(k)
    {
        this->marginalGainThreshold = ( deltaZero / ( 2 * k )) * pow(1 + epsilon, numBucket);
    }

    size_t getUtility() 
    {
        return localCovered.popcount();
    }

    int getTotalCovered()
    {
        return foundSeeds;
    }

    bool attemptInsert(const std::pair<int, std::unordered_set<int>>& element) 
    {
        std::unordered_set<int> temp;
        for (const int & e : element.second) {
            if (!localCovered.get(e))
                temp.insert(e);
        }
        if (temp.size() >= marginalGainThreshold) {
            this->seedSetAndMarginalGain[foundSeeds].first = element.first;
            this->seedSetAndMarginalGain[foundSeeds].second = temp.size();
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

    public:
    StreamingRandGreedIEngine(size_t theta, int deltaZero, int k, double epsilon, size_t totalBuckets) 
    {

        for (size_t i = 0; i < totalBuckets; i++) 
        {
            ThresholdBucket b(theta, deltaZero, k, epsilon, i);
            buckets.push_back(b);
        }

        this->theta = theta;
        this->k = k;
    }

    ~StreamingRandGreedIEngine() {}

    int insertElements(const std::vector<std::pair<int, std::unordered_set<int>>>& elements) 
    {
        int totalInsertions = 0;

        #pragma omp parallel for reduction(+:totalInsertions)
        for (size_t t = 0; t < buckets.size(); t++) 
        {
            for (const auto & element : elements) 
            {
                if (buckets[t].getTotalCovered() < this->k)
                {
                    totalInsertions += int(buckets[t].attemptInsert(element));
                }
            }
        }

        return totalInsertions;
    }

    template <typename A>
    std::pair<std::vector<unsigned int>, int> stream(CommunicationEngine<A> cEngine, int world_size) 
    {
        // While finished_data_streams < total_streams {
            // TODO: Iterate over your vector of receive buffers, std::vector<std::pair<status, int*>> 
            //  for every buffer that has data (test status) read it into another vector of type 
            //  std::vector<std::pair<int, std::unordered_map>> that holds your more permanent data. 

            // TODO: While reading stream data, if the tag sent with the data is a 1, this means stop 
            //  reading from this stream. 

            // TODO: Using the data read into permanent storage from the last step, in parallel iterate
            //  through all buckets in this.buckets. Loop again over all received data and try to instert
            //  each element into the bucket (which the bucket should independently handle). 

            // TODO: Buckets only insert IFF they contain less than k elements and the inserted element
            //  meets this bucket's threshold requirement. 

            // TODO: Handle stop conditions. For now the simplist (and least efficent) case will be assumed
            //  that the algorithm only stops when all streams have stopped. 
        // }

        // std::vector<std::pair<MPI_Request, int*>> receive_buffers(world_size);

        // int stoppingCount = 0; 
        // int streamedElements = 0;
        // while (stoppingCount <= (numBucket * k) || streamedElements <= (m * k)) {

        //     stoppingCount += insert(element);
        //     streamedElements++;
        // }

        // int max_marginal = 0;
        // int solution_bucket = 0;
        
        // for (auto & t : buckets) {

        //     if (max_marginal > t.getUtility()) {
        //         max_marginal = t.getUtility();
        //         solution_bucket++; 
        //     }
        // } 

        std::pair<std::vector<unsigned int>, int> result;

        // for (auto & s : (buckets.begin() + solution_bucket).seedSetAndMarginalGain)
        //     result.insert(std::make_pair(s.first, s.second));


        return result;
    }

};