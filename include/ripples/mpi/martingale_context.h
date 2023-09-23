#include "ripples/generate_rrr_sets.h"
#include "ripples/mpi/find_most_influential.h"
#include "ripples/utility.h"
#include "ripples/mpi/imm.h"
#include "ripples/imm_execution_record.h"
#include "ripples/imm.h"

template <
    typename GraphTy,
    typename ConfTy,
    typename RRRGeneratorTy,
    typename diff_model_tag,
    typename execution_tag
>
class MartingaleContext {
    private:
    DefaultSampler<GraphTy, diff_model_tag, RRRGeneratorTy, execution_tag> &sampler;
    const std::vector<ApproximatorContext> &approximators;
    OwnershipManager<GraphTy> &ownershipManager;

    const GraphTy &G;
    const ConfTy &CFG;
    ripples::IMMExecutionRecord &record;
    const CommunicationEngine<GraphTy> &cEngine;
    int optimalExecutionPath = -1;

    const double l;
    size_t previousTheta;

    std::map<int, std::vector<int>> solutionSpace;
    std::vector<size_t> old_sampling_sizes;
    ripples::RRRsetAllocator<typename GraphTy::vertex_type> allocator;
    TransposeRRRSets<GraphTy> tRRRSets;
    const std::vector<int> vertexToProcess;

    TimerAggregator &timeAggregator;

    ripples::RRRsetAllocator<typename GraphTy::vertex_type> GetAllocator() {
        #if defined ENABLE_MEMKIND
            return allocator(libmemkind::kinds::DAX_KMEM_PREFERRED);
        #elif defined ENABLE_METALL
            return metall_manager_instance().get_allocator();
        #else
            return ripples::RRRsetAllocator<typename GraphTy::vertex_type>();
        #endif
    }

    std::pair<std::vector<unsigned int>, int> runMartingaleRound(
        const size_t thetaPrime,
        const size_t previousTheta
    ) {
        // delta will always be (theta / 2) / world_size

        size_t delta = ((thetaPrime - previousTheta) / this->cEngine.GetSize() + 1);

        if ((thetaPrime - previousTheta) / this->cEngine.GetSize() < 0)
        {
            // std::cout << "DELTA ERROR" << std::endl;
            exit(1);
        }

        this->record.ThetaPrimeDeltas.push_back(delta);

        // std::cout << "before sampling, " << this->previousTheta << ", " << delta << std::endl;
        this->timeAggregator.samplingTimer.startTimer();
        this->sampler.addNewSamples(this->tRRRSets, this->previousTheta, delta);
        this->timeAggregator.samplingTimer.endTimer();    

        // std::cout << "before redistribution" << std::endl;
        this->timeAggregator.allToAllTimer.startTimer();  
        this->ownershipManager.redistributeSeedSets(this->tRRRSets, this->solutionSpace, delta);
        this->timeAggregator.allToAllTimer.endTimer();

        int kprime = int(CFG.alpha * (double)CFG.k);

        // std::cout << "before seed selection using kprime of " << kprime << std::endl;
        const std::map<int, std::vector<int>> &localSpace = this->solutionSpace;
        size_t localTheta = thetaPrime + this->cEngine.GetSize();

        if (this->approximators.size() == 1) {
            this->record.OptimalExecutionPath = 0;
            return this->approximators[0].getBestSeeds(localSpace, kprime, localTheta);
        }

        if (this->optimalExecutionPath != -1) {
            return this->approximators[this->optimalExecutionPath].getBestSeeds(localSpace, kprime, localTheta);
        }

        std::vector<std::future<std::pair<std::vector<unsigned int>, unsigned int>>> executionPaths;
        for (const auto & approximator : this->approximators) {
            executionPaths.push_back(
                std::async(
                    std::launch::async, 
                    [&approximator, &localSpace, &kprime, &localTheta] {
                        return approximator.getBestSeeds(localSpace, kprime, localTheta);
                    }
                )
            );
        }

        while (true) {
            for (size_t i = 0; i < executionPaths.size(); i++) {
                auto & f = executionPaths[i];
                std::future_status status;
                status = f.wait_for(std::chrono::seconds(0));
                if (status == std::future_status::ready) {
                    this->record.OptimalExecutionPath = i;
                    this->optimalExecutionPath = i;
                    return f.get();
                }
            }
        }
    }

    std::pair<std::vector<unsigned int>, unsigned int> runMartingaleRound(size_t theta) {
        if (this->previousTheta > theta)
        {
            std::cout << "invalid theta value, smaller than previous theta. Theta = " << theta << ", previous theta = " << this->previousTheta << std::endl;
            exit(1);
        }

        auto res = this->runMartingaleRound(
            theta, this->previousTheta
        );

        this->previousTheta = theta;

        return res;
    }

    inline double logBinomial(size_t n, size_t k) {
        return n * log(n) - k * log(k) - (n - k) * log(n - k);
    }

    inline size_t getTheta(
        double epsilon, double l, size_t k, double LB, size_t num_nodes
    ) {
        if (LB == 0) return 0;

        k = std::min(k, num_nodes/2);
        double term1 = 0.6321205588285577;  // 1 - 1/e
        double alpha = sqrt(l * std::log(num_nodes) + std::log(2));
        double beta = sqrt(term1 * (logBinomial(num_nodes, k) + l * std::log(num_nodes) + std::log(2)));
        double lamdaStar = 2 * num_nodes * (term1 * alpha + beta) * (term1 * alpha + beta) * pow(epsilon, -2);
        return lamdaStar / LB;
    }

    ssize_t getThetaPrime(
        ssize_t x, double epsilonPrime, double l, size_t k, size_t num_nodes 
    ) {
        k = std::min(k, num_nodes/2);
        return (2 + 2. / 3. * epsilonPrime) *
            (l * std::log(num_nodes) + logBinomial(num_nodes, k) +
            std::log(std::log2(num_nodes))) *
            std::pow(2.0, x) / (epsilonPrime * epsilonPrime);
    }

    static size_t DEBUG_countEmpty(const std::map<int, std::vector<int>> &localSolutionSpace){
	size_t emptyCount = 0;
	for (const auto & e : localSolutionSpace) {
	    if (e.second.size() == 0){
		emptyCount++;
	    }
	}

	return emptyCount;
    }

    public:
    MartingaleContext(
        DefaultSampler<GraphTy, diff_model_tag, RRRGeneratorTy, execution_tag> &sampler,
        OwnershipManager<GraphTy> &ownershipManager,
        const std::vector<ApproximatorContext> &approximators,

        const GraphTy &input_G, 
        const ConfTy &input_CFG,
        const double input_l,
        ripples::IMMExecutionRecord &record,
        const CommunicationEngine<GraphTy> &cEngine,
        TimerAggregator &timeAggregator
    ) 
        : 
            approximators(approximators), 
            ownershipManager(ownershipManager), 
            sampler(sampler),

            G(input_G),
            CFG(input_CFG), 
            tRRRSets(G.num_nodes()), 
            l(input_l * (1 + 1 / std::log2(G.num_nodes()))),
            vertexToProcess(cEngine.DistributeVertices(input_CFG.use_streaming, input_G)),
            old_sampling_sizes(input_G.num_nodes(), 0),
            record(record),
            cEngine(cEngine),
            timeAggregator(timeAggregator) {
        this->allocator = this->GetAllocator();
        this->previousTheta = 0;

        for (size_t i = 0; i < vertexToProcess.size(); i++)
        {
            if (vertexToProcess[i] == this->cEngine.GetRank())
            {
                this->solutionSpace.insert({i, std::vector<int>()});
            }
        }
    }

    std::vector<unsigned int> approximateInfMax() 
    {
        double LB = 0;
        double epsilonPrime = 1.4142135623730951 * this->CFG.epsilon;
        size_t thetaPrime = 0;

        std::pair<std::vector<unsigned int>, int> seeds;

        ////// MARTINGALE LOOP //////

        auto start = std::chrono::high_resolution_clock::now();

        for (ssize_t x = 1; x < std::log2(G.num_nodes()); ++x) 
        {
            // Equation 9
            thetaPrime = this->getThetaPrime(
                x, epsilonPrime, this->l, this->CFG.k, this->G.num_nodes()
            );

            // std::cout << "global theta: " << thetaPrime << std::endl;

            seeds = this->runMartingaleRound(thetaPrime);

            // f is the fraction of RRRsets covered by the seeds / the total number of RRRSets (in the current iteration of the martingale loop)
            // this has to be a global value, if one process succeeds and another fails it will get stuck in communication (the algorithm will fail). 
            double f;

            if (this->cEngine.GetRank() == 0) {
                f = (double)(seeds.second) / thetaPrime;
                // std::cout << "thetaprime: " << thetaPrime << std::endl;
            }
            
            // mpi_broadcast f(s)
            this->timeAggregator.broadcastTimer.startTimer();
            this->cEngine.distributeF(&f);
            this->timeAggregator.broadcastTimer.endTimer();

            // std::cout << "seeds.second: (covered RRRSet IDs) = " << seeds.second << " , thetaPrme: " << thetaPrime << " , f = " << f << std::endl;
            if (f >= std::pow(2, -x)) {
                // std::cout << "Fraction " << f << std::endl;
                LB = (G.num_nodes() * f) / (1 + epsilonPrime);
                // spdlog::get("console")->info("Lower Bound {}", LB);
                break;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();

        size_t theta = this->getTheta(this->CFG.epsilon, this->l, this->CFG.k, LB, this->G.num_nodes());

        this->record.ThetaEstimationTotal = end - start;

        this->record.Theta = theta;

        if (this->cEngine.GetRank() == 0)
        {
            // spdlog::get("console")->info("Previous ThetaPrime: {}, current Theta: {}", thetaPrime, theta);
        }

        std::pair<std::vector<typename GraphTy::vertex_type>, int> bestSeeds;

        if (thetaPrime >= theta) 
        {
            bestSeeds = seeds; 
        }
        else 
        {
            bestSeeds = this->runMartingaleRound(theta);
        }

        return bestSeeds.first;
    }
};
