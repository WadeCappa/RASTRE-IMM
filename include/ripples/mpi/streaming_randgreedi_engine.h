#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <iostream> 
#include <utility>
#include <chrono>
#include "bitmask.h"


class thresholdBuckets {
    
    private:
    ripples::Bitmask<int> localCovered;
    std::unordered_set<std::pair<int,int>> seedSetAndMarginalGain;
    double marginalGainThreshold;

    public:
    thresholdBuckets(size_t theta, int deltaZero, int k, double epsilon, size_t numBucket) {

        localCovered(ripples::Bitmask<int> covered(theta));
        marginalGainThreshold = ( deltaZero / ( 2 * k )) * ((1 + epsilon) ^ numBucket );

    }

    size_t getUtility() {
        return localCovered.popcount();
    }

};


template <typename GraphTy>
class StreamingRandGreedIEngine {


    private:
    std::vector<thresholdBuckets> buckets;

    public:
    StreamingRandGreedIEngine(size_t theta, int deltaZero, int k, double epsilon, size_t totalBuckets) {

        for (size_t i = 0; i < totalBuckets; i++) {

            thresholdBuckets b(theta, deltaZero, k, epsilon, i);
            buckets.push_back(b);

        }
        

    }
    ~StreamingRandGreedIEngine() {}

    int insert(std::pair<int, std::unordered_set<int>> element) {

        int totalInsertions = 0;

    #pragma omp parallel for reduction(+:totalInsertions)
        for (size_t t = 0; t < buckets.size(); t++) {

            if (buckets[t].seedSetAndMarginalGain.size() < k){

                std::unordered_set<int> temp;
                for (int e : element.second) {
                    if (!localCovered.get(e))
                        temp.insert(e)
                }
                if (temp.size() >= marginalGainThreshold) {
                    buckets[t].seedSetAndMarginalGain.insert(std::make_pair(element.first, temp.size()));                    
                    totalInsertions += 1;
                }

            }
            else continue;
    
        }

        return totalInsertions;

    }

    std::unordered_set<std::pair<int,int>> stream() {


        //TODO:: Read from Receive buffer into more permanent storage


        int stoppingCount = 0; 
        int streamedElements = 0;
        while (stoppingCount <= (numBucket * k) || streamedElements <= (m * k)) {

            stoppingCount += insert(element);
            streamedElements++;
        }

        int max_marginal = 0;
        int solution_bucket = 0;
        
        for (auto & t : buckets) {

            if (max_marginal > t.getUtility()) {
                max_marginal = t.getUtility();
                solution_bucket++; 
            }
        } 

        std::unordered_set<std::pair<int,int>> result;

        for (auto & s : (buckets.begin() + solution_bucket).seedSetAndMarginalGain)
            result.insert(std::make_pair(s.first, s.second));


        return result;
    }

}