#include "ripples/generate_rrr_sets.h"
#include "ripples/imm_execution_record.h"
#include "ripples/mpi/find_most_influential.h"
#include "ripples/utility.h"
#include "ripples/mpi/imm.h"
#include "ripples/imm_execution_record.h"
#include "ripples/imm.h"

template <
    typename GraphTy,
    typename ConfTy,
    typename RRRGeneratorTy
>
class MartingaleRunner {
    private:
    SamplerContext<GraphTy> &sampler;
    const ApproximatorContext &approximator;
    SeedSetAggregator<GraphTy> &aggregator;

    const GraphTy &G;
    const ConfTy &CFG;
    RRRGeneratorTy &gen;
    ripples::IMMExecutionRecord &record;
    const CommunicationEngine<GraphTy> &cEngine;

    const double l;
    size_t previousTheta;

    std::map<int, std::vector<int>> aggregateSets;
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

    std::pair<std::vector<unsigned int>, int> martingaleRound(
        const size_t thetaPrime,
        const size_t previousTheta
    ) {
        // delta will always be (theta / 2) / world_size

        std::pair<std::vector<unsigned int>, size_t> approximated_solution;

        auto timeRRRSets = ripples::measure<>::exec_time([&]() 
        {
            // size_t delta = localThetaPrime - this->RR_sets;
            // size_t delta = this->RR_sets == 0 ? thetaPrime / this->cEngine.GetSize() : (thetaPrime / 2) / this->cEngine.GetSize();
            size_t delta = ((thetaPrime - previousTheta) / this->cEngine.GetSize() + 1);

            if ((thetaPrime - previousTheta) / this->cEngine.GetSize() < 0)
            {
                std::cout << "DELTA ERROR" << std::endl;
                exit(1);
            }

            this->record.ThetaPrimeDeltas.push_back(delta);

            // spdlog::get("console")->info("sampling ...");

            this->timeAggregator.samplingTimer.startTimer();

            // std::cout << this->RR_sets << ", " << delta << ", " << thetaPrime << ", " << localThetaPrime << std::endl;
            this->sampler.addNewSamples(this->tRRRSets, this->previousTheta, delta);
            this->timeAggregator.samplingTimer.endTimer();    

            this->timeAggregator.allToAllTimer.startTimer();  

            // spdlog::get("console")->info("distributing samples with AllToAll ...");
            this->aggregator.aggregateSeedSets(this->tRRRSets, delta, this->aggregateSets);
            
            this->timeAggregator.allToAllTimer.endTimer();

            // spdlog::get("console")->info("seed selection ...");

            int kprime = int(CFG.alpha * (double)CFG.k);

            approximated_solution = this->approximator.getBestSeeds(kprime, thetaPrime + this->cEngine.GetSize());
        });
        
        this->record.ThetaEstimationGenerateRRR.push_back(timeRRRSets);
        auto timeMostInfluential = ripples::measure<>::exec_time([&]() { });
        this->record.ThetaEstimationMostInfluential.push_back(timeMostInfluential);

        return approximated_solution;
    }

    std::pair<std::vector<unsigned int>, unsigned int> runSeedSelection(size_t theta) {
        if (this->previousTheta > theta)
        {
            std::cout << "invalid theta value, smaller than previous theta. Theta = " << theta << ", previous theta = " << this->previousTheta << std::endl;
            exit(1);
        }

        auto res = this->martingaleRound(
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

    public:
    MartingaleRunner(
        SamplerContext<GraphTy> &sampler,
        SeedSetAggregator<GraphTy> &aggregator,
        const ApproximatorContext &approximator,

        const GraphTy &input_G, 
        const ConfTy &input_CFG,
        const double input_l,
        RRRGeneratorTy &gen,
        ripples::IMMExecutionRecord &record,
        const CommunicationEngine<GraphTy> &cEngine,
        TimerAggregator &timeAggregator
    ) 
        : 
            approximator(approximator), 
            aggregator(aggregator), 
            sampler(sampler),

            G(input_G),
            CFG(input_CFG), 
            tRRRSets(G.num_nodes()), 
            l(input_l * (1 + 1 / std::log2(G.num_nodes()))),
            vertexToProcess(cEngine.DistributeVertices(input_CFG.use_streaming, input_G)),
            old_sampling_sizes(input_G.num_nodes(), 0),
            gen(gen),
            record(record),
            cEngine(cEngine),
            timeAggregator(timeAggregator) {
        this->allocator = this->GetAllocator();
        this->previousTheta = 0;

        for (size_t i = 0; i < vertexToProcess.size(); i++)
        {
            if (vertexToProcess[i] == this->cEngine.GetRank())
            {
                this->aggregateSets.insert({i, std::vector<int>()});
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

            std::cout << "global theta: " << thetaPrime << std::endl;

            seeds = this->runSeedSelection(thetaPrime);

            // f is the fraction of RRRsets covered by the seeds / the total number of RRRSets (in the current iteration of the martigale loop)
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
            bestSeeds = this->runSeedSelection(theta);
        }

        // if (this->CFG.output_diagnostics == true)
        // {
        //     this->OutputDiagnosticData();
        // }

        return bestSeeds.first;
    }
};