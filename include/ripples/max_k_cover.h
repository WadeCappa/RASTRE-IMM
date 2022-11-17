#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <mutex>
#include <iostream>
#include "bitmask.h"
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

    std::vector<int> max_cover_lazy_greedy(std::unordered_map<int, std::unordered_set<int>> data, int k, int theta) 
    {
        CompareMaxHeap<int> cmp;

        std::vector<std::pair<int, std::unordered_set<int>>> data_vec	= std::vector<std::pair<int, std::unordered_set<int>>>(data.begin(), data.end());
        std::priority_queue<std::pair<int, std::unordered_set<int>>,
                            std::vector<std::pair<int, std::unordered_set<int>>>,
                            decltype(cmp)> pq(data_vec.begin(), data_vec.end());

        ripples::Bitmask<int> covered(theta);
        
        std::vector<int> result(k, -1);
        int count = 0;


        while(count < k) {
            
            auto l = pq.top();
            pq.pop();            
            
            // remove RR IDs from l that are already covered. 
            std::unordered_set<int> temp;
            for (int e: l.second) {
                if (covered.get(e)) {
                    temp.insert(e);
                }
            }
            
            
            for (int e: temp) {
                l.second.erase(e); 
            }
            
            
            //calculate marginal gain
            auto marginal_gain = l.second.size();

            //calculate utiluty of next best
            auto r = pq.top();
            
            
            // if marginal of l is better than r's utility, l is the current best     
            if (marginal_gain >= r.second.size()) {
                result[count] = l.first;
                count += 1;
                
                for (auto e : l.second) {
                    covered.set(e);
                }
            }
            // else push l's marginal into the heap 
            else {
                pq.push(l);
                
            }
        }

        return result;
    }
};