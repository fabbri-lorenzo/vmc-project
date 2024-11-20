#include "statistics.inl"

namespace vmcp {

// Helper function for accumulators
void FillAcc_(AccumulatorSet &dataAcc, const std::vector<Energy> &energyVector) {
    for (const auto &enValue : energyVector) {
        dataAcc(enValue.val);
    }
}

// Core helper method for Blocking
BlockingResult BlockingAnalysis_(std::vector<Energy> const &energies) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    // Check numEnergies is a power of 2 using bitwise AND between the number n and (n - 1)
    assert(numEnergies > 0 && (numEnergies & (numEnergies - 1)) == 0);
    assert(numEnergies > 1);
    std::vector<IntType> blockSizeList;
    std::vector<FPType> meansList;
    std::vector<FPType> stdDevsList;

    for (IntType blockSize = 1; blockSize <= numEnergies / 2; blockSize *= 2) {
        int numOfBlocks = (numEnergies + blockSize - 1) / blockSize;
        // Evaluate mean of each block
        std::vector<Energy> blockMeans = FillStatBlocks_<1>(energies, numEnergies, blockSize, numOfBlocks);
        // Evaluate second moment (mean squared) of each block
        std::vector<Energy> blockSecondM = FillStatBlocks_<3>(energies, numEnergies, blockSize, numOfBlocks);
        AccumulatorSet accBlocksMeans;
        AccumulatorSet accBlocksSecondM;
        FillAcc_(accBlocksMeans, blockMeans);
        FillAcc_(accBlocksSecondM, blockSecondM);
        // Statistics
        blockSizeList.push_back(blockSize);
        FPType mean = getStatValue_<1>(accBlocksMeans);
        meansList.push_back(mean);
        FPType stdDev =
            std::sqrt((getStatValue_<1>(accBlocksSecondM) - std::pow(mean, 2)) / (numOfBlocks - 1));
        stdDevsList.push_back(stdDev);
    }
    assert(blockSizeList.size() == meansList.size());
    assert(blockSizeList.size() == stdDevsList.size());
    return BlockingResult{blockSizeList, meansList, stdDevsList};
}

FPType BlockingOut(std::vector<Energy> const &energies) {
    BlockingResult blockingResult = BlockingAnalysis_(energies);
    // Find maximum value of standard dev.
    auto maxIt = std::max_element(blockingResult.stdDevs.begin(), blockingResult.stdDevs.end());
    FPType bestStdDev = *maxIt;
    return bestStdDev;
}

// Helper function for Bootstrapping
std::vector<Energy> GenerateBootstrapSample_(std::vector<Energy> const &energies,
                                             UIntType const &numEnergies_, RandomGenerator &gen,
                                             std::uniform_int_distribution<> &dist) {
    std::vector<Energy> sample;
    // Resample with replacement
    std::generate_n(std::back_inserter(sample), numEnergies_,
                    [&]() { return energies[static_cast<UIntType>(dist(gen))]; });
    return sample;
}

BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    std::uniform_int_distribution<> uIntDis(0, (numEnergies - 1));
    ConfInterval confInterval;
    // Calculate mean of each sample vector and place into bootstrapMeans
    std::vector<Energy> bootstrapMeans =
        FillBootstrapVec_<1>(energies, numEnergies, numSamples, gen, uIntDis);
    // Calculate second moment (squared mean) of each sample vector and place into bootstrapMeans
    std::vector<Energy> bootstrapSecondM =
        FillBootstrapVec_<3>(energies, numEnergies, numSamples, gen, uIntDis);
    AccumulatorSet accBlocksMeans;
    AccumulatorSet accBlocksSecondM;
    FillAcc_(accBlocksMeans, bootstrapMeans);
    FillAcc_(accBlocksSecondM, bootstrapSecondM);
    // Calculate statistics of each means
    FPType meanOfMeans = getStatValue_<1>(accBlocksMeans);
    FPType stdDev =
        std::sqrt((getStatValue_<1>(accBlocksSecondM) - std::pow(meanOfMeans, 2)) / (numEnergies - 1));
    confInterval.min = meanOfMeans - stdDev * zScore;
    confInterval.max = meanOfMeans + stdDev * zScore;
    return BootstrapResult{meanOfMeans, stdDev, confInterval};
}
} // namespace vmcp
