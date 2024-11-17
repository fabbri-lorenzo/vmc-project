#include "statistics.inl"

namespace vmcp {

BlockingResult BlockingAnalysis(std::vector<Energy> const &energies) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    std::vector<IntType> blockSizeList;
    std::vector<FPType> meansList;
    std::vector<FPType> stdDevsList;

    for (IntType blockSize = 1; blockSize <= numEnergies / 2; blockSize += 2) {
        int numOfBlocks = (numEnergies + blockSize - 1) / blockSize;
        std::vector<FPType> blockMeans(static_cast<UIntType>(numOfBlocks));
        std::vector<FPType> blockSecondM(static_cast<UIntType>(numOfBlocks));
        // Evaluate mean of each block
        FillStatBlocks_<1>(blockMeans, energies, numEnergies, blockSize);
        // Evaluate second moment (mean squared) of each block
        FillStatBlocks_<3>(blockSecondM, energies, numEnergies, blockSize);
        AccumulatorSet accBlocksMeans;
        AccumulatorSet accBlocksSecondM;
        FillAcc_<FPType>(accBlocksMeans, blockMeans);
        FillAcc_<FPType>(accBlocksSecondM, blockSecondM);
        // Statistics
        blockSizeList.push_back(blockSize);
        FPType mean = getStatValue_<1>(accBlocksMeans);
        meansList.push_back(mean);
        FPType stdDev =
            std::sqrt((getStatValue_<1>(accBlocksSecondM) - std::pow(mean, 2)) / (numOfBlocks - 1));
        stdDevsList.push_back(stdDev);

        if (blockSize == 1) {
            blockSize = 0;
        }
    }
    assert(blockSizeList.size() == meansList.size());
    assert(blockSizeList.size() == stdDevsList.size());
    return BlockingResult{blockSizeList, meansList, stdDevsList};
}

// Helper function for Bootstrapping
std::vector<Energy> GenerateBootstrapSample_(std::vector<Energy> const &energies,
                                             UIntType const &numEnergies_, RandomGenerator &gen,
                                             std::uniform_int_distribution<> &uIntDis) {
    std::vector<Energy> sample(numEnergies_);
    // Resample with replacement
    std::generate(sample.begin(), sample.end(),
                  [&]() { return energies[static_cast<UIntType>(uIntDis(gen))]; });
    return sample;
}

BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    std::uniform_int_distribution<> uIntDis(0, (numEnergies - 1));
    ConfInterval confInterval;
    std::vector<FPType> bootstrapMeans(numSamples);
    std::vector<FPType> bootstrapSecondM(numSamples);
    // Calculate mean of each sample vector and place into bootstrapMeans
    FillBootstrapVec_<1>(bootstrapMeans, energies, numEnergies, gen, uIntDis);
    // Calculate second moment (squared mean) of each sample vector and place into bootstrapMeans
    FillBootstrapVec_<3>(bootstrapSecondM, energies, numEnergies, gen, uIntDis);
    AccumulatorSet accBlocksMeans;
    AccumulatorSet accBlocksSecondM;
    FillAcc_<FPType>(accBlocksMeans, bootstrapMeans);
    FillAcc_<FPType>(accBlocksSecondM, bootstrapSecondM);
    // Calculate statistics of each means
    FPType meanOfMeans = getStatValue_<1>(accBlocksMeans);
    FPType stdDev =
        std::sqrt((getStatValue_<1>(accBlocksSecondM) - std::pow(meanOfMeans, 2)) / (numEnergies - 1));
    confInterval.min = meanOfMeans - stdDev * zScore;
    confInterval.max = meanOfMeans + stdDev * zScore;
    return BootstrapResult{meanOfMeans, stdDev, confInterval};
}
} // namespace vmcp
