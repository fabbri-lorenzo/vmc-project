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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <ranges>
#include <vector>

namespace vmcp {

//  Critical value for confidence level of 95% in Gaussian distributed data
FPType constexpr zScore = 1.96;

// Helper function to get the desired statistic value from the accumulator
template <int Stat>
FPType getStatValue_(AccumulatorSet &dataAcc) {
    if constexpr (Stat == STAT_MEAN) {
        return boost::accumulators::mean(dataAcc);
    } else if constexpr (Stat == STAT_VARIANCE) {
        return boost::accumulators::variance(dataAcc);
    } else if constexpr (Stat == STAT_SECONDMOMENT) {
        return boost::accumulators::moment<2>(dataAcc);
    } else {
        static_assert(Stat == STAT_MEAN || Stat == STAT_VARIANCE || Stat == STAT_SECONDMOMENT,
                      "Unsupported statistic.");
    }
}

// Helper function for accumulators
template <typename T>
void FillAcc_(AccumulatorSet &dataAcc, const std::vector<T> &energyVector) {
    if constexpr (std::is_same_v<T, Energy>) {
        for (const auto &enValue : energyVector) {
            dataAcc(enValue.val);
        }
    } else if constexpr (std::is_same_v<T, FPType>) {
        for (const auto &enValue : energyVector) {
            dataAcc(enValue);
        }
    } else {
        static_assert(std::is_same_v<T, Energy> || std::is_same_v<T, FPType>, "Unsupported type.");
    }
}

// Helper function for Blocking, fills blocks-vectors with desired statistic
template <int Stat>
void FillStatBlocks_(std::vector<FPType> &blockStats, std::vector<Energy> const &energies,
                     IntType numEnergies, IntType blockSize) {
    AccumulatorSet accBlock;
    int rem = numEnergies % blockSize;
    IntType currentBlock = 0;
    std::generate(blockStats.begin(), blockStats.end(),
                  [&energies, &numEnergies, &accBlock, &blockSize, &currentBlock, &rem]() mutable {
                      accBlock = {}; // Reset accumulator
                      auto start = energies.begin() + currentBlock * blockSize;
                      auto end = (currentBlock * blockSize + blockSize <= numEnergies) ? start + blockSize
                                                                                       : start + rem;
                      FillAcc_<Energy>(accBlock, std::vector<Energy>(start, end));
                      ++currentBlock;
                      return getStatValue_<Stat>(accBlock);
                  });
}

// Helper function for Bootstrapping
std::vector<Energy> GenerateBootstrapSample_(std::vector<Energy> const &energies,
                                             UIntType const &numEnergies_, RandomGenerator &gen,
                                             std::uniform_int_distribution<> &uIntDis);

// Helper function for Bootstrapping, fills sample-vectors with desired statistic
template <int Stat>
void FillBootstrapVec_(std::vector<FPType> &bootstrapVector, std::vector<Energy> const &energies,
                       int const &numEnergies, RandomGenerator &gen,
                       std::uniform_int_distribution<> &uIntDis) {
    AccumulatorSet accSample;
    std::generate(bootstrapVector.begin(), bootstrapVector.end(),
                  [&accSample, &energies, &numEnergies, &gen, &uIntDis]() {
                      accSample = {};
                      std::vector<Energy> bootstrapSample = GenerateBootstrapSample_(
                          energies, static_cast<UIntType>(numEnergies), gen, uIntDis);
                      FillAcc_<Energy>(accSample, bootstrapSample);
                      return getStatValue_<Stat>(accSample);
                  });
}
} // namespace vmcp

#endif