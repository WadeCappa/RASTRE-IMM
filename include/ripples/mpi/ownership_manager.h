

template <typename GraphTy>
class OwnershipManager {
    private:
    public:
    virtual void redistributeSeedSets(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        std::map<int, std::vector<int>> &aggregateSets,
        const size_t delta
    ) = 0;
};

template <typename GraphTy>
class AllToAllOwnershipManager : public OwnershipManager<GraphTy> {
    private: 
        const CommunicationEngine<GraphTy> &cEngine;
        const std::vector<int> vertexToProcess;
        std::vector<size_t> oldSamplingSizes;

    public:

    AllToAllOwnershipManager(
        size_t num_nodes,
        const CommunicationEngine<GraphTy> &cEngine,
        const std::vector<int> vertexToProcess
    ) : cEngine(cEngine), vertexToProcess(vertexToProcess), oldSamplingSizes(num_nodes, 0) {

    }
    
    void redistributeSeedSets(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        std::map<int, std::vector<int>> &aggregateSets,
        const size_t delta
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