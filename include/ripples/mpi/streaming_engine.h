#include <vector>
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>
#include <mutex>
#include <iostream> 
#include <utility>
#include <chrono>
#include <math.h>
#include <algorithm>
#include <deque>
#include "ripples/bitmask.h"

#include <thread>
#include <future>

#include <cmath>

typedef struct origin {
    int source;
    int seed;
} Origin;

typedef struct candidateSet {
    unsigned int vertex;
    std::vector<unsigned int> covered;
} CandidateSet;

class ThresholdBucket 
{
    private:
    ripples::Bitmask<int> localCovered;
    std::vector<std::pair<int,int>> seeds;
    double marginalGainThreshold;
    int k;
    ssize_t theta;

    public:
    ThresholdBucket(ssize_t theta, int deltaZero, int k, double epsilon, size_t numBucket) 
        : localCovered(theta)
    {
        this->theta = theta;
        this->marginalGainThreshold = ( (double)deltaZero / (double)( 2 * k )) * (double)std::pow(1 + epsilon, numBucket);
        this->k = k;
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
        return this->seeds.size();
    }

    bool attemptInsert(const CandidateSet& element) 
    {
        if (this->seeds.size() >= this->k)
        {
            return false;
        }
       
        std::vector<int> temp;
        for (const int e : element.covered) {
            if (!localCovered.get(e))
                temp.push_back(e);
        }

        if (temp.size() >= this->marginalGainThreshold) {
            for (const int e : temp) {
                localCovered.set(e);
            }

            this->seeds.push_back(std::make_pair(element.vertex, temp.size()));             
            return true;
        }

        return false;
    }
};

class BucketController 
{
    private:
    std::vector<ThresholdBucket*>* buckets = new std::vector<ThresholdBucket*>();
    int k;
    ssize_t theta;
    double epsilon;
    
    int deltaZero = -1;    
    int fullBuckets = 0;

    public:
    void CreateBuckets(size_t deltaZero)
    {
        this->deltaZero = deltaZero;

        int num_buckets = (int)(0.5 + [](double val, double base) {
            return log2(val) / log2(base);
        }((double)k, (1+this->epsilon)));

        std::cout << "deltazero: " << this->deltaZero << ", number of buckets: " << num_buckets << std::endl;

        for (int i = 0; i < num_buckets + 1; i++)
        {
            this->buckets->push_back(new ThresholdBucket(this->theta, this->deltaZero, this->k, this->epsilon, i));
        }
    }

    size_t GetNumberOfBuckets()
    {
        return this->buckets->size();
    }

    std::vector<ThresholdBucket*>* GetBuckets()
    {
        return this->buckets;
    }

    bool Initialized()
    {
        return this->buckets->size() > 0;
    }

    bool AllBucketsFull()
    {
        // return this->buckets.size() == 0 ? false : this->fullBuckets == this->buckets.size();
        return false;
    }

    std::pair<std::vector<unsigned int>, int> GetBestSeeds()
    {
        size_t max_covered = 0;
        int max_covered_index = 0;

        for (int i = 0; i < this->buckets->size(); i++)
        {
            size_t bucket_utility = this->buckets->at(i)->getUtility();
            if (bucket_utility > max_covered) {
                max_covered = bucket_utility;
                max_covered_index = i;
            }
        }
        
        std::vector<unsigned int> seeds;
        for (const auto p : this->buckets->at(max_covered_index)->getSeeds()) {
            seeds.push_back(p.first);
        }

        std::cout << "selected bucket " << max_covered_index << " with size of " << this->buckets->at(max_covered_index)->getSeeds().size() << std::endl;

        return std::make_pair(seeds, max_covered);
    }

    BucketController(int k, ssize_t theta, double epsilon)
    {
        this->theta = theta;
        this->k = k;
        this->epsilon = epsilon;
    }

    ~BucketController() 
    {
        for (auto b : *(this->buckets))
        {
            delete b;
        }

        delete this->buckets;
    }
};

template <typename GraphTy>
class StreamingRandGreedIEngine 
{
    private: 
    int* buffer;
    MPI_Request request;
    BucketController buckets;
    const CommunicationEngine<GraphTy> &cEngine;

    int k;
    int kprime;
    int active_senders;
    int first_values_from_senders;
    int number_of_senders;
    size_t theta;
    double epsilon;
    TimerAggregator& timer;

    std::vector<std::vector<CandidateSet>> element_matrix;
    std::vector<std::vector<unsigned int>> local_seed_sets;
    std::vector<std::pair<int, Origin>> availability_index;
    std::vector<unsigned int> local_utilities;

    static void GetSeedSet(int* data, std::vector<unsigned int>& seed_set)
    {
        int negatives = 0;
        int* e = data;
        for (; negatives != 2; e++)
        {
            if (negatives == 1 && *e != -1)
            {
                seed_set.push_back(*e);
            }

            if (*e == -1)
            {
                negatives++;
            }
        }
    }

    static CandidateSet ExtractElement(int* data)
    {
        std::vector<unsigned int> received_data;

        for (int* e = data + 1; *(e) != -1; e++) 
        {
            received_data.push_back(*e);
        }

        return (CandidateSet){(unsigned int)*data, received_data};
    }

    void ResetBuffer(int* buffer)
    {
        this->cEngine.QueueReceive(buffer, this->theta, this->request);
    }

    void HandleStatus(const MPI_Status& status)
    {
        if (status.MPI_TAG == 0)
        {
            this->first_values_from_senders++;
        }

        if (status.MPI_TAG > this->kprime - 1)
        {
            this->active_senders--;
        }
    }

    int* ReceiveNextSend(MPI_Status& status)
    {
        this->timer.receiveTimer.startTimer();
        MPI_Wait(&(this->request), &status);
        this->timer.receiveTimer.endTimer();

        return this->buffer;
    }

    public:
    StreamingRandGreedIEngine(
        int k, int kprime, size_t theta, 
        double epsilon, int number_of_senders, 
        const CommunicationEngine<GraphTy> &cEngine,
        TimerAggregator& timer
    ) 
        : buckets(k, theta, epsilon), cEngine(cEngine), timer(timer)
    {
        this->active_senders = number_of_senders;
        this->k = k;
        this->kprime = kprime;
        this->epsilon = epsilon;
        this->number_of_senders = number_of_senders;
        this->first_values_from_senders = 0;
        this->theta = theta;

        std::cout << "theta = " << theta << std::endl;

        this->buffer = new int[theta];
    
        this->element_matrix = std::vector<std::vector<CandidateSet>>(
            this->number_of_senders, 
            std::vector<CandidateSet>(
                this->kprime,
                (CandidateSet){0, std::vector<unsigned int>()}
            )
        );

        this->local_seed_sets = std::vector<std::vector<unsigned int>>(this->number_of_senders);

        this->availability_index = std::vector<std::pair<int, Origin>>(
            this->number_of_senders * this->kprime, 
            std::make_pair(0, (Origin){-1,-1})
        );

        this->local_utilities = std::vector<unsigned int>(this->number_of_senders);
    }

    ~StreamingRandGreedIEngine()
    {

    }

    std::pair<std::vector<unsigned int>, int> Stream()
    {
        int buckets_initialized = 0;
        int dummy_value = 0;
        int kill_processors = 0;

        int threads = omp_get_max_threads();
        omp_set_nested(2);

        this->timer.totalGlobalStreamTimer.startTimer();

        # pragma omp parallel num_threads(2)
        {
            if (omp_get_thread_num() == 0) // receiver
            {
                size_t maxVal = 0;

                // std::cout << "RECEIVING WITH THREAD ID " << omp_get_thread_num() << std::endl;

                for (int i = 0; i < (this->number_of_senders * this->kprime); i++)
                {                   
                    MPI_Status status;
                    this->ResetBuffer(this->buffer);
                    int* data = this->ReceiveNextSend(status);

                    this->timer.processingReceiveTimer.startTimer();

                    int tag = status.MPI_TAG > this->kprime - 1 ? this->kprime - 1 : status.MPI_TAG;
                    int source = status.MPI_SOURCE - 1;

                    // std::cout << "tag: " << tag << ", source: " << source << std::endl;

                    if (status.MPI_TAG > this->kprime - 1)
                    {
                        // std::cout << "got last seed" << std::endl;

                        local_utilities[source] = status.MPI_TAG;
                        this->GetSeedSet(data, local_seed_sets[source]);
                    }

                    element_matrix[source][tag] = this->ExtractElement(data);

                    this->timer.processingReceiveTimer.endTimer();

                    // std::cout << "extracted elements" << std::endl;

                    maxVal = std::max(maxVal, element_matrix[source][tag].covered.size());

                    availability_index[i].second.source = source;
                    availability_index[i].second.seed = tag;

                    this->timer.atomicUpdateTimer.startTimer();

                    #pragma omp atomic 
                    availability_index[i].first++;
                    
                    this->timer.atomicUpdateTimer.endTimer();

                    this->HandleStatus(status);

                    if (this->first_values_from_senders == this->number_of_senders && buckets_initialized == 0)
                    {
                        this->timer.initBucketTimer.startTimer();
                        this->buckets.CreateBuckets(maxVal);
                        this->timer.initBucketTimer.endTimer();

                        #pragma omp atomic 
                        buckets_initialized++;
                    }

                    // std::cout << "handled status" << std::endl;
                }

                std::cout << "killing processors, waiting for them to exit..." << std::endl;
            }
            else // processor
            {   
                while (true) 
                {
                    # pragma omp atomic
                    dummy_value++;
                    
                    if (buckets_initialized == 1)
                    {
                        break;
                    }
                }

                std::cout << "exited waiting loop..." << std::endl;

                // build mapping of each thread to number of buckets
                auto buckets = this->buckets.GetBuckets();
                int buckets_per_thread = std::ceil((double)buckets->size() / ((double)threads - (double)1));
                std::cout << "number of buckets; " << buckets->size() << ", threads; " << threads-1 << ", buckets per thread; " << buckets_per_thread << std::endl;
                std::vector<std::vector<ThresholdBucket*>> bucketMap(threads-1);

                int bucket_index = 0;
                for (int i = 0; i < bucketMap.size(); i++ )
                {
                    int buckets_added = 0;
                    while (bucket_index < buckets->size() && buckets_added < buckets_per_thread)
                    {
                        bucketMap[i].push_back(buckets->at(bucket_index));
                        bucket_index++;
                        buckets_added++;
                    }
                }

                std::cout << "starting to process elements..." << std::endl;

                this->timer.max_k_globalTimer.startTimer();

                #pragma omp parallel for num_threads(threads-1)
                for (int i = 0; i < bucketMap.size(); i++)
                {
                    const std::vector<ThresholdBucket*>& thread_buckets = bucketMap[i];
                    int local_dummy_value = 0;

                    for (int local_received_index = 0; local_received_index < availability_index.size(); local_received_index++)
                    {
                        while (true) 
                        {
                            if (availability_index[local_received_index].first == 1)
                            {
                                break;
                            }

                            # pragma omp atomic
                            local_dummy_value++;
                        }

                        // std::cout << "inserting seed " << local_received_index << std::endl;

                        for (const auto & b : thread_buckets)
                        {
                            auto p = availability_index[local_received_index].second;
                            b->attemptInsert(element_matrix[p.source][p.seed]);
                        }
                    }
                }
                
                this->timer.max_k_globalTimer.endTimer();
            }
        }

        this->timer.totalGlobalStreamTimer.endTimer();

        auto bestSeeds = this->buckets.GetBestSeeds();

        for (int i = 0; i < local_utilities.size(); i++)
        {
            // std::cout << "streaming utility is " << bestSeeds.second << std::endl;
            // std::cout << "looking at " << local_utilities[i] << " from local process " << i+1 << std::endl;
            if (local_utilities[i] > bestSeeds.second)
            {
                bestSeeds.second = local_utilities[i];
                bestSeeds.first = local_seed_sets[i];
            }
        }

        std::cout << "number of seeds: " << bestSeeds.first.size() << std::endl;

        return bestSeeds;
    }
};