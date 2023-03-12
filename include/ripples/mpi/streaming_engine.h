#include <vector>
#include <unordered_set>
#include <unordered_map>
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

using ElementList = std::vector<std::pair<int, std::unordered_set<int>*>>;

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

    bool attemptInsert(const std::pair<int, std::unordered_set<int>*>& element) 
    {
        if (this->seeds.size() == this->k)
        {
            return false;
        }
       
        std::unordered_set<int> temp;
        for (const int e : *(element.second)) {
            if (e > this->theta)
                std::cout << "ELEMENT " << e << " BIGGER THAN THETA; " << this->theta << std::endl;
            else if (!localCovered.get(e))
                temp.insert(e);
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
    std::vector<ThresholdBucket*> buckets;
    int k;
    ssize_t theta;
    double epsilon;
    
    int deltaZero = -1;    
    int fullBuckets = 0;

    public:
    void CreateBuckets(ElementList* current_elements)
    {
        // calculate deltaZero
        size_t maxval = 0;
        for (const auto & r : *current_elements)
        {
            maxval = std::max(maxval, r.second->size());
        }

        this->deltaZero = maxval;

        int num_buckets = (int)(0.5 + [](double val, double base) {
            return log2(val) / log2(base);
        }((double)k, (1+this->epsilon)));

        std::cout << "deltazero: " << this->deltaZero << ", number of buckets: " << num_buckets << std::endl;

        for (int i = 0; i < num_buckets + 1; i++)
        {
            this->buckets.push_back(new ThresholdBucket(this->theta, this->deltaZero, this->k, this->epsilon, i));
        }
    }

    bool Initialized()
    {
        return this->buckets.size() > 0;
    }

    void ProcessData(const ElementList* elements, const size_t threads)
    {
        // std::cout << "processing " << elements->size() << " elements..." << std::endl;

        #pragma omp parallel for num_threads(threads)
        for (size_t t = 0; t < buckets.size(); t++)
        {
            for (int i = 0; i < elements->size(); i++) 
            {
                buckets[t]->attemptInsert(elements->at(i));
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

        for (int i = 0; i < this->buckets.size(); i++)
        {
            size_t bucket_utility = this->buckets[i]->getUtility();
            if (bucket_utility > max_covered) {
                max_covered = bucket_utility;
                max_covered_index = i;
            }
        }
        
        std::vector<unsigned int>* seeds = new std::vector<unsigned int>();
        for (const auto p : this->buckets[max_covered_index]->getSeeds()) {
            seeds->push_back(p.first);
        }

        return std::make_pair(*seeds, max_covered);
    }

    BucketController(int k, ssize_t theta, double epsilon)
    {
        this->theta = theta;
        this->k = k;
        this->epsilon = epsilon;
    }

    ~BucketController() 
    {
        for (const auto b : this->buckets)
        {
            delete b;
        }
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

    static std::pair<int, std::unordered_set<int>*> ExtractElement(int* data)
    {
        std::unordered_set<int>* received_data = new std::unordered_set<int>();

        for (int* e = data + 1; *(e) != -1; e++) 
        {
            received_data->insert(*e);
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

    bool HandleStatus(MPI_Status& status)
    {
        bool buckets_initialized = false;

        if (status.MPI_TAG == 0)
        {
            this->first_values_from_senders++;

            if (this->first_values_from_senders == this->world_size)
            {
                this->buckets.CreateBuckets(this->elements);
                buckets_initialized = true;
            }
        }

        if (status.MPI_TAG == this->k - 1)
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
        this->active_senders = world_size;
        this->k = k;
        this->epsilon = epsilon;
        this->world_size = world_size;
        this->first_values_from_senders = 0;
        this->theta = theta;

        this->elements = new ElementList();

        this->buffer = new int[theta];
        this->ResetBuffer();
    }

    std::pair<std::vector<unsigned int>, int> Stream(Timer* timer)
    {
        MPI_Status status;
        std::mutex* lock = new std::mutex(); 
        bool streaming_finished = false;
        bool buckets_initialized = false;

        size_t threads = omp_get_num_threads();

        # pragma omp parallel num_threads(2) shared(lock, buckets_initialized)
        {
            if (omp_get_thread_num() == 0) // receiver
            {
                for (int i = 0; i < (this->world_size * this->k) && this->active_senders > 0; i++)
                {
                    MPI_Wait(this->request, &status);

                    lock->lock();

                    buckets_initialized = this->HandleStatus(status);

                    this->elements->push_back(this->ExtractElement(this->buffer));

                    lock->unlock();

                    if (i != this->world_size * this->k - 1)
                    {
                        this->ResetBuffer();
                    }
                }

                streaming_finished = true;
            }
            else // processor
            {   
                ElementList* local_elements = 0;
                int number = 0;

                while (true) 
                {
                    lock->lock();

                    if (this->buckets.Initialized())
                    {
                        lock->unlock();
                        break;
                    }

                    lock->unlock();
                }

                std::cout << "exited while loop, no longer waiting" << std::endl;

                while (streaming_finished == false)
                {
                    lock->lock();

                    if (this->elements->size() > 0)
                    {
                        local_elements = this->elements;
                        this->elements = new ElementList(); 
                    }

                    lock->unlock();

                    if (local_elements != 0)
                    {
                        timer->startTimer();
                        this->buckets.ProcessData(local_elements, threads - 1);
                        timer->endTimer();
                        
                        delete local_elements;
                        local_elements = 0;
                    }
                }
            }
        }

        return this->buckets.GetBestSeeds();
    }
};