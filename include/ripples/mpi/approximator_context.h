

class ApproximatorContext {
    private:
    public:

    ApproximatorContext() {}

    std::pair<std::vector<unsigned int>, unsigned int> getBestSeeds(
        unsigned int kprime,
        size_t theta
    ) const {
        return std::make_pair(std::vector<unsigned int>{0,1,2,3,4,5,6,7,8,9}, 30000);
    }
};