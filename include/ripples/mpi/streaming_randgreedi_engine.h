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


class ThresholdBucket 
{
    private:
    ripples::Bitmask<int> localCovered;
    std::vector<std::pair<int,int>> seeds;
    double marginalGainThreshold;
    int foundSeeds;
    int k;

    public:
    ThresholdBucket(size_t theta, int deltaZero, int k, double epsilon, size_t numBucket) 
        : localCovered(theta), seeds(k)
    {
        this->marginalGainThreshold = ( deltaZero / ( 2 * k )) * std::pow(1 + epsilon, numBucket);
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
        return foundSeeds;
    }

    bool attemptInsert(const std::pair<int, std::unordered_set<int>*>& element) 
    {
        std::unordered_set<int> temp;
        for (const int & e : *(element.second)) {
            if (!localCovered.get(e))
                temp.insert(e);
        }
        if (temp.size() >= this->marginalGainThreshold && this->foundSeeds < this->k) {
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
    typedef struct receiveObject {
        int* data;
        MPI_Request request;
        int receive_count;
    } ReceiveObject;

    std::vector<ThresholdBucket*> buckets;
    size_t theta;
    int k;
    int world_size;
    int deltaZero = -1;
    double epsilon;
    std::vector<ReceiveObject> receive_buffers;
    std::deque<std::pair<int, std::unordered_set<int>*>> elements;


    void createBuckets()
    {
        // calculate deltaZero
        size_t maxval = 0;
        for (const auto & r : this->elements)
        {
            maxval = std::max(maxval, r.second->size());
        }

        this->deltaZero = maxval;

        int num_buckets = (int)(0.5 + [](double base, double val) {
            return log2(val) / log2(base);
        }((1+this->epsilon), ((k-1)*(this->deltaZero/(2*this->k)))));

        for (int i = 0; i < num_buckets + 1; i++)
        {
            this->buckets.push_back(new ThresholdBucket(this->theta, this->deltaZero, this->k, this->epsilon, i));
        }
    }

    int handleNewElements()
    {
        int received_elements = 0;
        int local_process_sent_first_seed = 0;

        for (int i = 0; i < this->world_size; i++) 
        {
            local_process_sent_first_seed += this->receive_buffers[i].receive_count > 0 ? 1 : 0;
            
            int flag;
            MPI_Status status;
            MPI_Test(&(this->receive_buffers[i].request), &flag, &status);
            if (flag == 1) {
                received_elements++;

                std::unordered_set<int>* received_data = new std::unordered_set<int>();

                for (int e = *(this->receive_buffers[i].data + 1); e != -1; e = *(&e + 1)) 
                {
                    received_data->insert(e);
                }

                this->elements.push_front(std::make_pair(*(this->receive_buffers[i].data), received_data));

                MPI_Irecv (
                    this->receive_buffers[i].data,
                    this->theta, MPI_INT, i, 0, 
                    MPI_COMM_WORLD, 
                    &(this->receive_buffers[i].request)
                );
            }
        }

        if (local_process_sent_first_seed == this->k && this->deltaZero == -1)
        {
            this->createBuckets();
        }

        return received_elements;
    }

    void processData()
    {
        if (this->deltaZero == -1)
            return;

        while (this->elements.size() > 0)
        {
            std::pair<int, std::unordered_set<int>*> next_element = this->elements[0];
            this->elements.pop_back();

            #pragma omp parallel for
            for (size_t t = 0; t < buckets.size(); t++) 
            {        
                if (buckets[t]->getTotalCovered() < this->k)
                {
                    int(buckets[t]->attemptInsert(next_element));
                }                
            }

            delete next_element.second;
        }
    }

    bool allBucketsFull()
    {
        return false;
    }

    public:
    StreamingRandGreedIEngine(int k, int theta, double epsilon, int world_size)
        : receive_buffers(world_size)
    {
        this->theta = theta;
        this->k = k;
        this->epsilon = epsilon;
        this->world_size = world_size;

        // initialization step
        for (int i = 0; i < world_size; i++) 
        {
            // In the future this could be a bitmap to save space and lower communciation times
            this->receive_buffers[i].data = new int[this->theta];

            MPI_Irecv (
                this->receive_buffers[i].data,
                this->theta, MPI_INT, i, 0,
                MPI_COMM_WORLD,
                &(this->receive_buffers[i].request)
            );
        }
    }

    ~StreamingRandGreedIEngine() {}

    std::pair<std::vector<unsigned int>, int> stream() 
    {
        int received_elements = 0;
        while (received_elements < (this->k * this->world_size) && !this->allBucketsFull()) 
        {
            received_elements += this->handleNewElements();
            this->processData();
        }

        std::pair<std::vector<unsigned int>, int> result;

        int max_covered = 0;
        int max_covered_index = 0;

        for (int i = 0; i < this->buckets.size(); i++)
        {
            if (max_covered > this->buckets[i]->getUtility()) {
                max_covered = this->buckets[i]->getUtility();
                max_covered_index = i;
            }
        }
        
        std::vector<unsigned int>* seeds = new std::vector<unsigned int>();
        for (const auto & p : this->buckets[max_covered_index]->getSeeds()) {
            seeds->push_back(p.first);
        }

        return std::make_pair(*seeds, max_covered);
    }
};