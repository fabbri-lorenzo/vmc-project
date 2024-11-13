//
//
// Contains the definition of the templates declared in statistics.hpp
// This file is supposed to be #included at the end of statistics.hpp and nowhere else
// It is just a way to improve the readability of statistics.hpp
//

// The functions that have a trailing underscore would, in a regular project, be (declared and) defined in a
// .cpp file, therefore making them unreachable for the user
// Since they are templated, they must be defined in a header

#ifndef VMCPROJECT_STATISTICS_INL
#define VMCPROJECT_STATISTICS_INL

#include "statistics.hpp"

#include <boost/accumulators/accumulators.hpp>
// LF TODO: this file should not be in framework dir usually, check why it got installed there
#include <algorithm>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics.hpp>
#include <cmath>
#include <fstream>
#include <numeric>
#include <vector>

namespace vmcp {
// Boost library type
typedef boost::accumulators::accumulator_set<
    FPType, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance,
                                          boost::accumulators::tag::min, boost::accumulators::tag::max>>
    AccumulatorSet;

BlockingResult BlockingAnalysis(std::vector<Energy> const &energies) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    std::vector<IntType> blockSizeList;
    std::vector<FPType> meansList;
    std::vector<FPType> stdDevsList;

    for (IntType blockSize = 1; blockSize <= numEnergies / 2; blockSize += 2) {
        // LF TODO: this is actually the number of FULL blocks -> incomplete name?
        int numOfBlocks = numEnergies / blockSize;
        std::vector<FPType> blockMeans(static_cast<UIntType>(numOfBlocks));
        IntType block = 0;

        std::transform(blockMeans.begin(), blockMeans.end(), blockMeans.begin(),
                       [&energies, blockSize, block](FPType &blockMean) mutable {
                           blockMean =
                               std::accumulate(
                                   energies.begin() + static_cast<DiffType>(block * blockSize),
                                   energies.begin() + static_cast<DiffType>((block + 1) * blockSize),
                                   Energy{0}, [](Energy e1, Energy e2) { return Energy{e1.val + e2.val}; })
                                   .val /
                               blockSize;
                           ++block;
                           return blockMean;
                       });
        // Cases in which numEnergies is not a multiple of blockSize
        int reminder = numEnergies % blockSize;
        if (reminder != 0) {
            blockMeans.push_back(
                std::accumulate(energies.begin() + static_cast<DiffType>(numOfBlocks),
                                energies.begin() + static_cast<DiffType>(numOfBlocks + reminder), Energy{0},
                                [](Energy e1, Energy e2) { return Energy{e1.val + e2.val}; })
                    .val /
                reminder);
        }

        // Creating the accumulator set with mean and variance features
        AccumulatorSet accBlocks;
        for (const auto &value : blockMeans) {
            accBlocks(value);
        }
        blockSizeList.push_back(blockSize);
        meansList.push_back(boost::accumulators::mean(accBlocks));
        stdDevsList.push_back(std::sqrt(boost::accumulators::variance(accBlocks)));

        if (blockSize == 1) {
            blockSize = 0;
        }
    }
    assert(blockSizeList.size() == meansList.size());
    assert(blockSizeList.size() == stdDevsList.size());
    return BlockingResult{blockSizeList, meansList, stdDevsList};
}

// Helper function for bootstrapping
template <typename AccumulatorFunc>
FPType GenerateBootstrapSample(std::vector<Energy> const &energies, UIntType const &numEnergies_,
                               RandomGenerator &gen, std::uniform_int_distribution<> &uIntDis,
                               AccumulatorFunc &accumulator, std::string statistic) {
    std::vector<FPType> sample(numEnergies_);
    // Resample with replacement
    std::generate(sample.begin(), sample.end(),
                  [&]() { return energies[static_cast<UIntType>(uIntDis(gen))].val; });

    // Accumulate the values in the sample
    for (const auto &enValue : sample) {
        accumulator(enValue);
    }

    if (statistic == "mean") {
        return boost::accumulators::mean(accumulator);
    }
    if (statistic == "variance") {
        return boost::accumulators::variance(accumulator);
    } else {
        throw std::invalid_argument("Please provide a correct statistics label");
    }
}

BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen) {
    int const numEnergies = static_cast<IntType>(std::ssize(energies));
    assert(numEnergies > 1);
    ConfInterval confInterval;
    std::vector<FPType> bootstrapMeans(numSamples);
    std::vector<FPType> bootstrapVars(numSamples);
    std::uniform_int_distribution<> uIntDis(0, (numEnergies - 1));

    // Define accumulators for means and variances
    AccumulatorSet accBlocksMeans;
    AccumulatorSet accBlocksVars;

    // Fill bootstrap means and variances using the helper
    std::generate(bootstrapMeans.begin(), bootstrapMeans.end(), [&]() {
        return GenerateBootstrapSample(energies, static_cast<UIntType>(numEnergies), gen, uIntDis,
                                       accBlocksMeans, "mean");
    });

    std::generate(bootstrapVars.begin(), bootstrapVars.end(), [&]() {
        return GenerateBootstrapSample(energies, static_cast<UIntType>(numEnergies), gen, uIntDis,
                                       accBlocksVars, "variance");
    });

    // Clear accumulators and recalculate overall means and vars
    AccumulatorSet accBlocksMeansCleared;
    for (const auto &value : bootstrapMeans) {
        accBlocksMeansCleared(value);
    }
    AccumulatorSet accBlocksVarsCleared;
    for (const auto &value : bootstrapVars) {
        accBlocksVarsCleared(value);
    }

    FPType meanOfMeans = boost::accumulators::mean(accBlocksMeansCleared);
    FPType stdDevOfMeans = std::sqrt(boost::accumulators::mean(accBlocksVarsCleared) / (numEnergies - 1));

    confInterval.min = boost::accumulators::min(accBlocksMeansCleared);
    confInterval.max = boost::accumulators::max(accBlocksMeansCleared);

    return BootstrapResult{meanOfMeans, stdDevOfMeans, confInterval};
}
} // namespace vmcp

#endif
