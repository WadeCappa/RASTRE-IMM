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

using ElementList = std::vector<std::pair<int, std::vector<int>*>>;

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

    bool attemptInsert(const std::pair<int, std::vector<int>*>& element) 
    {
        if (this->seeds.size() == this->k)
        {
            return false;
        }
       
        std::vector<int> temp;
        for (const int e : *(element.second)) {
            if (!localCovered.get(e))
                temp.push_back(e);
        }

        if (temp.size() >= this->marginalGainThreshold) {
            for (const int e : temp) {
                localCovered.set(e);
            }

            this->seeds.push_back(std::make_pair(element.first, temp.size()));             
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

    void ProcessData(const ElementList* elements)
    {
        // std::cout << "processing " << elements->size() << " elements..." << std::endl;

        #pragma omp parallel for 
        for (size_t t = 0; t < buckets->size(); t++)
        {
            for (int i = 0; i < elements->size(); i++) 
            {
                buckets->at(t)->attemptInsert(elements->at(i));
            }
        }
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

class StreamingRandGreedIEngine 
{
    private: 
    int* buffer;
    MPI_Request* request;
    BucketController buckets;

    int k;
    int active_senders;
    int first_values_from_senders;
    int world_size;
    ssize_t theta;
    double epsilon;

    ElementList* elements;

    static std::pair<int, std::vector<int>*> ExtractElement(int* data)
    {
        std::vector<int>* received_data = new std::vector<int>();

        for (int* e = data + 1; *(e) != -1; e++) 
        {
            received_data->push_back(*e);
        }

        return std::make_pair(*data, received_data);
    }

    void ResetBuffer()
    {
        MPI_Irecv(
            this->buffer,
            this->theta,
            MPI_INT,
            MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            MPI_COMM_WORLD,
            this->request
        );
    }

    bool HandleStatus(MPI_Status& status, int kprime)
    {
        bool buckets_initialized = false;
        if (status.MPI_TAG == 0)
        {
            this->first_values_from_senders++;

            if (this->first_values_from_senders == this->world_size)
            {
                buckets_initialized = true;
            }
        }

        if (status.MPI_TAG >= kprime - 1)
        {
            // this means that last element from a process has been sent
            this->active_senders--;
        }

        return buckets_initialized;
    }

    public:
    StreamingRandGreedIEngine(int k, ssize_t theta, double epsilon, int world_size) 
        : buckets(k, theta, epsilon)
    {
        this->active_senders = world_size - 1;
        this->k = k;
        this->epsilon = epsilon;
        this->world_size = world_size;
        this->first_values_from_senders = 0;
        this->theta = theta;

        this->elements = new ElementList();

        this->buffer = new int[theta];

        this->request = new MPI_Request();
        
        this->ResetBuffer();
    }

    ~StreamingRandGreedIEngine()
    {
        delete this->elements;
        delete this->buffer;
        delete this->request;
    }

    std::pair<std::vector<unsigned int>, int> Stream(TimerAggregator* timer, int kprime)
    {
        MPI_Status status;
        int buckets_initialized = 0;
        int dummy_value = 0;
        int kill_processors = 0;
        bool streaming_finished = false;

        int threads = omp_get_max_threads();

        std::vector<std::vector<std::pair<int, std::vector<int>*>>> element_matrix(
            this->world_size, 
            std::vector<std::pair<int, std::vector<int>*>>(
                kprime,
                std::make_pair(-1, (std::vector<int>*)0)
            )
        );

        std::vector<std::pair<int, std::pair<int,int>>> availability_index(
            this->world_size * kprime, 
            std::make_pair(0, std::make_pair(-1,-1))
        );

        std::vector<int> local_utilities(this->world_size);

        omp_set_nested(2);

        timer->totalGlobalStreamTimer.startTimer();

        # pragma omp parallel num_threads(2) shared(availability_index, element_matrix, buckets_initialized, dummy_value, kill_processors)
        {
            if (omp_get_thread_num() == 0) // receiver
            {
                size_t maxVal = 0;

                // std::cout << "RECEIVING WITH THREAD ID " << omp_get_thread_num() << std::endl;

                for (int i = 0; i < (this->world_size * kprime); i++)
                {
                    // TODO: Add all stop conidtions

                    // TODO: Create a method for communicating with sending processes
                    //  after all buckets have filled up.
                    
                    timer->receiveTimer.startTimer();
                    MPI_Wait(this->request, &status);
                    timer->receiveTimer.endTimer();

                    timer->processingReceiveTimer.startTimer();

                    int tag = status.MPI_TAG > this-> k - 1 ? kprime - 1 : status.MPI_TAG;
                    int source = status.MPI_SOURCE - 1;

                    if (status.MPI_TAG > kprime - 1)
                    {
                        local_utilities[source] = status.MPI_TAG;
                        // std::cout << "received utility of " << status.MPI_TAG << " from " << source << std::endl;
                    }

                    auto new_element = this->ExtractElement(this->buffer);

                    element_matrix[source][tag] = new_element;

                    maxVal = std::max(maxVal, new_element.second->size());

                    if (this->HandleStatus(status, kprime))
                    {
                        timer->initBucketTimer.startTimer();
                        this->buckets.CreateBuckets(maxVal);
                        timer->initBucketTimer.endTimer();

                        #pragma omp atomic 
                        buckets_initialized++;
                    }

                    availability_index[i].second.first = source;
                    availability_index[i].second.second = tag;

                    timer->processingReceiveTimer.endTimer();

                    timer->atomicUpdateTimer.startTimer();

                    #pragma omp atomic 
                    availability_index[i].first++;
                    
                    timer->atomicUpdateTimer.endTimer();


                    if (i != this->world_size * kprime - 1)
                    {
                        this->ResetBuffer();
                    }
                }

                // std::cout << "killing processors, waiting for them to exit..." << std::endl;

                streaming_finished = true;
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

                // std::cout << "exited waiting loop..." << std::endl;

                // build mapping of each thread to number of buckets
                auto buckets = this->buckets.GetBuckets();
                int buckets_per_thread = std::ceil((double)buckets->size() / ((double)threads - (double)1));
                std::cout << "number of buckets; " << buckets->size() << ", threads; " << threads-1 << ", buckets per thread; " << buckets_per_thread << std::endl;
                std::vector<std::vector<ThresholdBucket*>> bucketMap(threads-1, std::vector<ThresholdBucket*>());

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

                // std::cout << "starting to process elements..." << std::endl;

                timer->max_k_globalTimer.startTimer();

                #pragma omp parallel for num_threads(threads-1)
                for (int i = 0; i < bucketMap.size(); i++)
                {
                    auto thread_buckets = bucketMap[i];
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

                        for (auto & b : thread_buckets)
                        {
                            auto p = availability_index[local_received_index].second;
                            b->attemptInsert(element_matrix[p.first][p.second]);
                        }
                    }
                }
                timer->max_k_globalTimer.endTimer();
            }
        }

        timer->totalGlobalStreamTimer.endTimer();

        auto bestSeeds = this->buckets.GetBestSeeds();

        for (int i = 0; i < local_utilities.size(); i++)
        {
            // std::cout << "streaming utility is " << bestSeeds.second << std::endl;
            // std::cout << "looking at " << local_utilities[i] << " from local process " << i+1 << std::endl;
            if (local_utilities[i] > bestSeeds.second)
            {
                bestSeeds.second = local_utilities[i];
                std::vector<unsigned int> new_best_seeds;
                for (const auto & s : element_matrix[i])
                {
                    new_best_seeds.push_back(s.first);
                }

                bestSeeds.first = new_best_seeds;
            }
        }

        // std::cout << "number of seeds: " << bestSeeds.first.size() << std::endl;

        for (auto & r : element_matrix)
        {
            for (auto & c : r)
            {
                if (c.second != 0)
                {
                    delete c.second;
                }
            }
        }

        return bestSeeds;
    }
};