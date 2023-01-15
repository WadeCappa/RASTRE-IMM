#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <mutex>
#include <iostream>
#include "bitmask.h"
#include <utility>
#include <queue>

class MaxKCoverEngine 
{
    public:

    template <typename idTy>
    struct CompareMaxHeap {

        bool operator()(const std::pair<idTy, std::unordered_set<idTy>> a,
                        const std::pair<idTy, std::unordered_set<idTy>> b) {
            return a.second.size() < b.second.size();
        }
    };

    std::pair<std::vector<unsigned int>, int> max_cover(std::unordered_map<int, std::unordered_set<int>>& data, int k, int theta) {

        ripples::Bitmask<int> covered(theta);
        std::vector<unsigned int> result(k, -1);
        int totalUniqueCovered = 0;

        for (size_t i = 0; i < k; i++) {

            int max = 0;
            int max_key = -1;

            for( const auto& t : data ) {

                if (t.second.size() > max) {
                    max = t.second.size();
                    max_key = t.first;
                }

            }
            result[i] = max_key;

            for (int e: data[max_key]) {
                if (!covered.get(e)) {
                    covered.set(e);
                    totalUniqueCovered++;
                }
            }

            for( const auto& t : data ) {

                std::set<int> temp;
                if (t.first != max_key) {
                    for (int e: t.second) {
                        if (covered.get(e)) {
                            temp.insert(e);
                        }
                    }
                    for (int e: temp) {
                        data[t.first].erase(e); 
                    }

                }

            }

            data.erase(max_key);
        }
        return std::make_pair(result, totalUniqueCovered);
    }

    std::pair<std::unordered_set<unsigned int>, int> max_cover_lazy_greedy(std::unordered_map<int, std::unordered_set<int>>& data, int k, int theta) 
    {
        CompareMaxHeap<int> cmp;
        int totalRRRSetsUsed = 0;

        std::vector<std::pair<int, std::unordered_set<int>>> data_vec = std::vector<std::pair<int, std::unordered_set<int>>>(data.begin(), data.end());
        std::priority_queue<std::pair<int, std::unordered_set<int>>,
                            std::vector<std::pair<int, std::unordered_set<int>>>,
                            decltype(cmp)> pq(data_vec.begin(), data_vec.end());

        ripples::Bitmask<int> covered(theta);
        
        std::unordered_set<unsigned int> result;

        std::cout << "theta = " << theta << " , k = " << k << " , data.size() = " << data.size() << std::endl;
        
        while(result.size() < k && pq.size() > 1) {                        
            auto l = pq.top();
            pq.pop();       

            std::unordered_set<int> temp;

            // remove RR IDs from l that are already covered. 
            for (int e: l.second) {
                if (e > theta || e < 0) {
                    std::cout << "e is greater than theta, e = " << e << " , theta = " << theta << std::endl;
                }
                else if (covered.get(e)) {
                    temp.insert(e);
                }
            }            

            for (const int e : temp) 
            {
                l.second.erase(e); 
            }
            
            // calculate marginal gain
            auto marginal_gain = l.second.size();

            // calculate utiluty of next best
            auto r = pq.top(); 
            
            // if marginal of l is better than r's utility, l is the current best     
            if (marginal_gain >= r.second.size()) {
                result.insert(l.first);
                
                for (int e : l.second) {
                    if (e > theta || e < 0) {
                        std::cout << "e is greater than theta, e = " << e << " , theta = " << theta << std::endl;
                    }
                    else if (!covered.get(e)) {
                        totalRRRSetsUsed++;
                        covered.set(e);
                    }
                }
            }
            // else push l's marginal into the heap 
            else {
                pq.push(l);                
            }
        }

        return std::make_pair(result, totalRRRSetsUsed);
    }
};