

template <typename GraphTy>
class SeedSetAggregator {
    private:
    public:
    virtual void aggregateSeedSets(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        size_t delta, 
        std::map<int, std::vector<int>> &aggregateSets
    ) = 0;
};

template <typename GraphTy>
class AllToAllSeedAggregator : public SeedSetAggregator<GraphTy> {
    private: 
        const CommunicationEngine<GraphTy> &cEngine;
        const std::vector<int> vertexToProcess;
        std::vector<size_t> oldSamplingSizes;

    public:

    AllToAllSeedAggregator(
        size_t num_nodes,
        const CommunicationEngine<GraphTy> &cEngine,
        const std::vector<int> vertexToProcess
    ) : cEngine(cEngine), vertexToProcess(vertexToProcess), oldSamplingSizes(num_nodes, 0) {

    }
    
    void aggregateSeedSets(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        const size_t delta, 
        std::map<int, std::vector<int>> &aggregateSets
    ) override {
        this->cEngine.AggregateThroughAllToAll(
            tRRRSets,
            this->vertexToProcess,
            this->oldSamplingSizes,
            delta,
            aggregateSets
        );
    }
};