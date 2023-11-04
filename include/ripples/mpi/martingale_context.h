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

    static size_t choose(size_t n, size_t k) { // TODO: This logic will have rounding errors, see what the opimc repo did for this.
        size_t res = 1;
        for (size_t i = 0; i < k; i++) {
            res *= (n - (k - i)) / i; 
        }
    
        return res;
    }

    static size_t calculateThetaMax(size_t num_nodes, double approx, double eps, size_t k, double delta) { // TODO: this logic will have rounding errors, only convert to size_t at the very end.
        const double denominator = std::pow(eps, 2) * k;
        const double natural_log = std::log(6.0 / delta);
        const size_t n_choose_k = choose(num_nodes, k);
        const double numerator = natural_log * approx + std::sqrt(approx * (natural_log + std::log(n_choose_k)));
        return (2 * num_nodes * std::pow(numerator, 2)) / denominator;
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

    std::vector<unsigned int> useImm() 
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


    std::vector<unsigned int> useOpimc() 
    {
        const double epsilonPrime = 1.4142135623730951 * this->CFG.epsilon;
        const double delta = 1.0 / (double)this->G.num_nodes();
        const double approximation_guarantee = 1.0 - (1.0 / (double)std::exp(1.0));

        const double error_delta = (double)1 / (double)G.num_nodes();
        const size_t theta_max = calculateThetaMax(G.num_nodes(), approximation_guarantee, this->CFG.epsilon, this->CFG.k, delta);
        const size_t theta_0 = std::ceil(((double)theta_max * (double)std::pow(this->CFG.epsilon, 2) * (double)this->CFG.k) / (double)this->G.num_nodes());
        size_t global_theta = theta_0;
        
        // will be modified
        size_t local_theta = std::floor((double)theta_0 / (double)this->cEngine.GetSize());
        size_t local_theta_delta = local_theta;

        TransposeRRRSets<GraphTy> R1(this->G.num_nodes());
        TransposeRRRSets<GraphTy> R2(this->G.num_nodes());

        const int i_max = std::ceil(std::log2(theta_max / theta_0));

        this->timeAggregator.samplingTimer.startTimer();
        this->sampler.addNewSamples(R1, 0, local_theta_delta);
        this->sampler.addNewSamples(R2, 0, local_theta_delta);
        this->timeAggregator.samplingTimer.endTimer();

        for (int i = 0; i < i_max; i++) {
            record.ThetaPrimeDeltas.push_back(local_theta_delta);

            this->timeAggregator.allToAllTimer.startTimer();  
            this->ownershipManager.redistributeSeedSets(R1, this->solutionSpace, local_theta);
            this->timeAggregator.allToAllTimer.endTimer();  

            const int kprime = int(this->CFG.alpha * (double)(this->CFG.k));

            const auto s_star = this->approximators[0].getBestSeeds(this->solutionSpace, kprime, global_theta);
            
            const std::vector<unsigned int> seeds = this->cEngine.distributeSeedSet(s_star.first, kprime);
            const unsigned int local_R2_influence = R2.calculateInfluence(seeds); /// Evaluate the influence spread of a seed set on current generated RR sets
            const unsigned int global_R2_influence = this->cEngine.sumCoverage(local_R2_influence);

            double alpha;
            if (this->cEngine.GetRank() == 0) {
                const double a1 = std::log(2.0 / delta);
                const double a2 = std::log(2.0 / delta);

                const unsigned int infSelf = s_star.second; /// influence over R1
                const auto degVldt = (double)(global_R2_influence * global_theta) / (double)(this->G.num_nodes());

                const auto upperBound = (double)infSelf / (double)approximation_guarantee;

                const auto upperDegOPT = (double)(upperBound * global_theta) / (double)(this->G.num_nodes());
                const double sigma_super_l = std::pow(std::sqrt(degVldt + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0), 2) - a1 / 18.0;
                const double sigma_super_u = std::pow(std::sqrt(upperDegOPT + a2 / 2.0) + sqrt(a2 / 2.0), 2);
                alpha = sigma_super_l / sigma_super_u;
            }
            
            this->cEngine.distributeF(&alpha);

            if (alpha >= (approximation_guarantee - this->CFG.epsilon)) { // divide approx by 2 here? Check in with group.
                return s_star.first;
            }

            local_theta_delta = local_theta;
            this->timeAggregator.samplingTimer.startTimer();
            this->sampler.addNewSamples(R1, global_theta, local_theta_delta);
            this->sampler.addNewSamples(R2, global_theta, local_theta_delta);
            this->timeAggregator.samplingTimer.endTimer();
            local_theta = local_theta * 2;
            global_theta = global_theta * 2;
        }
    
        std::cout << "Error, exited for loop without returning, should not be possible" << std::endl;
        return std::vector<unsigned int>();
    }
};
