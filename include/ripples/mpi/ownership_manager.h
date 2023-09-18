

template <typename GraphTy>
class OwnershipManager {
    private: 
        const CommunicationEngine<GraphTy> &cEngine;
        const std::vector<int> &vertexToProcess;
        std::vector<size_t> oldSamplingSizes;

    public:

    OwnershipManager(
        size_t num_nodes,
        const CommunicationEngine<GraphTy> &cEngine,
        const std::vector<int> &vertexToProcess
    ) : cEngine(cEngine), vertexToProcess(vertexToProcess), oldSamplingSizes(num_nodes, 0) {
	int rank = cEngine.GetRank();
	size_t num_responsible = 0;
	for (const auto & e : vertexToProcess) {
	    if (e == rank) {
		num_responsible++;
	    }
	}
	std::cout << "rank " << rank << " is responsible for " << num_responsible << " vertices" << std::endl;
    }
    
    void redistributeSeedSets(
        const TransposeRRRSets<GraphTy> &tRRRSets, 
        std::map<int, std::vector<int>> &aggregateSets,
        const size_t delta
    ) {
        this->cEngine.AggregateThroughAllToAll(
            tRRRSets,
            this->vertexToProcess,
            this->oldSamplingSizes,
            delta,
            aggregateSets
        );
    }
};
