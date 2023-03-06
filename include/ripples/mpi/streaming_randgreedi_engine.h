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

class ThresholdBucket 
{
    private:
    ripples::Bitmask<int> localCovered;
    std::vector<std::pair<int,int>> seeds;
    double marginalGainThreshold;
    int k;

    public:
    ThresholdBucket(size_t theta, int deltaZero, int k, double epsilon, size_t numBucket) 
        : localCovered(theta)
    {
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
        std::unordered_set<int> temp;
        for (const int e : *(element.second)) {
            if (!localCovered.get(e))
                temp.insert(e);
        }
        if (temp.size() >= this->marginalGainThreshold && this->seeds.size() < this->k) {
            for (const int e : temp) {
                localCovered.set(e);
            }

            this->seeds.push_back(std::make_pair(element.first, temp.size()));             
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
        bool active;
    } ReceiveObject;

    std::vector<ThresholdBucket*> buckets;
    size_t theta;
    int k;
    int world_size;
    int deltaZero = -1;
    double epsilon;
    std::vector<std::vector<ReceiveObject>> receive_buffers;

    int active_senders;
    int fullBuckets = 0;

    void createBuckets(std::vector<std::pair<int, std::unordered_set<int>*>>& current_elements)
    {
        // calculate deltaZero
        size_t maxval = 0;
        for (const auto & r : current_elements)
        {
            maxval = std::max(maxval, r.second->size());
        }

        this->deltaZero = maxval;

        int num_buckets = (int)(0.5 + [](double base, double val) {
            return log2(val) / log2(base);
        }((1+this->epsilon), ((double)(k-1)*((double)this->deltaZero/(double)(2*k)))));

        for (int i = 0; i < num_buckets + 1; i++)
        {
            this->buckets.push_back(new ThresholdBucket(this->theta, this->deltaZero, this->k, this->epsilon, i));
        }
    }

    void processData(std::vector<std::pair<int, std::unordered_set<int>*>>& elements)
    {
        if (this->deltaZero == -1)
            return;

        for (auto & element: elements)
        {

            #pragma omp parallel for
            for (size_t t = 0; t < buckets.size(); t++) 
            {
                if (buckets[t]->getTotalCovered() < this->k)
                {
                    int(buckets[t]->attemptInsert(element));
                }                
            }

            // when you refactor this function to process the whole queue in parallel, delete this in the parent function, not here
            delete element.second;
        }
    }

    bool allBucketsFull()
    {
        // return this->buckets.size() == 0 ? false : this->fullBuckets == this->buckets.size();
        return false;
    }

    std::pair<std::vector<unsigned int>, int> getBestSeeds()
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

    public:
    StreamingRandGreedIEngine(int k, int theta, double epsilon, int world_size)
        : receive_buffers(world_size)
    {
        this->theta = theta;
        this->k = k;
        this->epsilon = epsilon;
        this->world_size = world_size;
        this->active_senders = world_size;

        // initialization step
        for (int i = 0; i < this->world_size; i++) 
        {
            this->receive_buffers[i] = std::vector<ReceiveObject>(this->k);

            for (int j = 0; j < this->k; j++)
            {
                // In the future this could be a bitmap to save space and lower communciation times
                this->receive_buffers[i][j].data = new int[this->theta];
                this->receive_buffers[i][j].active = true;

                MPI_Irecv (
                    this->receive_buffers[i][j].data,
                    this->theta, MPI_INT, i, j,
                    MPI_COMM_WORLD,
                    &(this->receive_buffers[i][j].request)
                );
            }
        }
    }

    static int handleNewElements(
        std::vector<std::pair<int, std::unordered_set<int>*>>& elements,
        std::vector<std::vector<ReceiveObject>>* receive_buffers,
        int& active_senders,
        const size_t theta,
        int world_size
    )
    {
        int received_elements = 0;
        int local_process_sent_first_seed = 0;

        for (int i = 0; i < world_size; i++) 
        {            
            for (int j = 0; j < receive_buffers->at(j).size(); j++)
            {
                if (receive_buffers->at(i)[j].active)
                {
                    int flag;
                    MPI_Status status;
                    MPI_Test(&(receive_buffers->at(i)[j].request), &flag, &status);

                    if (flag == 1) 
                    {
                        received_elements++;

                        std::unordered_set<int>* received_data = new std::unordered_set<int>();

                        for (int* e = receive_buffers->at(i)[j].data + 1; *(e) != -1; e++) 
                        {
                            received_data->insert(*e);
                        }

                        elements.push_back(std::make_pair(*(receive_buffers->at(i)[j].data), received_data));
                        delete receive_buffers->at(i)[j].data;
                        receive_buffers->at(i)[j].active = false;
                    }
                }
            }
        }

        return received_elements;
    }

    ~StreamingRandGreedIEngine() {}

    std::pair<std::vector<unsigned int>, int> stream(Timer* timer) 
    {
        int received_elements = 0;
        std::vector<std::pair<int, std::unordered_set<int>*>>* q = new std::vector<std::pair<int, std::unordered_set<int>*>>();

        while 
        (
            received_elements < (this->k * this->world_size)  
            && !this->allBucketsFull()  
            && this->active_senders > 0
        ) 
        {
            if (this->buckets.size() == 0)
            {
                int new_element_count = StreamingRandGreedIEngine::handleNewElements(
                    *q,
                    &(this->receive_buffers), 
                    this->active_senders,
                    this->theta,
                    this->world_size
                );

                received_elements += new_element_count;

                int received_buffers = 0;
                for (const auto & b : receive_buffers)
                {
                    received_buffers += b[0].active ? 0 : 1;
                }
                
                if (received_buffers == this->world_size && this->buckets.size() == 0)
                {
                    this->createBuckets(*q);
                }
            }
            else 
            {
                auto new_q = new std::vector<std::pair<int, std::unordered_set<int>*>>();

                std::vector<std::vector<ReceiveObject>>* buffers = &(this->receive_buffers);
                int senders = this->active_senders;
                int t = this->theta;
                int m = this->world_size;

                auto ret = std::async(
                    [&new_q, &buffers, &senders, &t, &m](){
                        return StreamingRandGreedIEngine::handleNewElements(
                            *new_q,
                            buffers, 
                            senders,
                            t,
                            m
                        );
                    } 
                );

                timer->startTimer();
                this->processData(*q);
                timer->endTimer();

                int new_element_count = ret.get();

                this->active_senders = senders;

                delete q;
                q = new_q;
                received_elements += new_element_count;
            }
        }

        return this->getBestSeeds();
    }
};