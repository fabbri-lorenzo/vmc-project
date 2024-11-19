#include "statistics.inl"

// FP: No include guard?

// FP: I like "average" more than "mean", mostly because it follows the already present convention
// Which means: should we nchange the convention and always use "mean"?

namespace vmcp {

// FP: Ask for the energies to be of size "power of 2"
BlockingResult BlockingAnalysis(std::vector<Energy> const &energies) {
    // FP: IntType
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    // FP: rename to BlockSizes, etc.
    std::vector<IntType> blockSizeList;
    std::vector<FPType> meansList;
    std::vector<FPType> stdDevsList;

    // FP: Shouldn't you do blockSize *= 2, like "Error estimates on averages of correlated data" says?
    // FP: Does the blocking analysis end here?
    // Isn't the end goal to obtain a better estimate for the error by looking at the "plateau" of the
    // variance OF THE AVERAGE?
    for (IntType blockSize = 1; blockSize <= numEnergies / 2; blockSize += 2) {
        // FP: Hey, this is an integer division
        int numOfBlocks = (numEnergies + blockSize - 1) / blockSize;
        // FP: Is the static_cast really necessary here?
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

// FP: uIntDis could be renamed to "dist", since it is clear that it is a distribution of unsigned integers
// Helper function for Bootstrapping
std::vector<Energy> GenerateBootstrapSample_(std::vector<Energy> const &energies,
                                             UIntType const &numEnergies_, RandomGenerator &gen,
                                             std::uniform_int_distribution<> &uIntDis) {
    std::vector<Energy> sample(numEnergies_);
    // Resample with replacement
    std::generate(sample.begin(), sample.end(),
                  [&]() { return energies[static_cast<UIntType>(uIntDis(gen))]; });
    // FP: You could avoid specifiying the size of the vector and use generate_n with std::back_inserter
    return sample;
}

// FP: Not a fan of numSamples being UInt, but probably you are following some convention of mine which must
// therefore be changed
BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen) {
    // FP: IntType
    // FP: Also you could avoid the static_cast, I believe, as assignment through = does not check overflow
    // etc
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    // FP: Is <> really necessary?
    // Maybe it is and I should just shut up
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
    // FP: You do not need to return the confidence interval, that can easily be calculated from the other
    // objects that you return
    // Also not returning it means that the user can use his preferred number of std devs for the conf int
    return BootstrapResult{meanOfMeans, stdDev, confInterval};
}

} // namespace vmcp

// FP: I would prefer that BootstrapAnalysis and Blocking Analysis returned a VMCResult object, since their
// purpose is to better estimate the error on the vmc energy

// FP: Also, I believe that using boost to do "just" averages and variances is kinda overkill
// The optimal algorithm for the variance is the usual one (implemented in AvgAndVar_)
// I know that it's two-pass, but I looked for better algorithms for estimating the average and found nothing,
// just one-pass algorithms that are less precise

// FP: Also I believe that we should add tests for one particle in 2D, and for 2 particles in 1D, just to see
// that nothing breaks
