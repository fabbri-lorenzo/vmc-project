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
    // FP: careful with these if else
    // Go look "dangling else"
    // Or maybe this is a perfectly safe syntax and I should study more and talk less
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
// FP: Descritpion does not say what the function does
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
// FP: Since the functions expects blockStats to be empty, you do not care about the previous content of
// blockStats Then is it really necessary to have it as input parameter? Why not return it from the function?
// FP: Really not a fan of the fact that all blocks have the same size, except (maybe) for the last one
// Could you ask for a vector of energies of opportune size as input?
// Just tear down the program (i.e. fail an assert) if you do not get what you want, like a 3 year old
template <int Stat>
void FillStatBlocks_(std::vector<FPType> &blockStats, std::vector<Energy> const &energies,
                     IntType numEnergies, IntType blockSize) {
    AccumulatorSet accBlock;
    // FP: IntType
    int rem = numEnergies % blockSize;
    IntType currentBlock = 0;
    // FP: This function is so insanely dangerous (not that you are insane, the function is)
    // It relies on a delicate balance between the size of blockStatsm numEnergies, and blockSize
    // I know that you only call this when you are sure that the size of blockStats is ...
    // But this constraint is not at all obvious is general, it requires someone to look inside the machinery
    // of your function to see how he should behave, which exactly the opposite of what functions are for
    std::generate(blockStats.begin(), blockStats.end(),
                  [&energies, &numEnergies, &accBlock, &blockSize, &currentBlock, &rem]() mutable {
                      // FP: Is mutable really necessary here?
                      // Also, I would just define accBlock inside this {}, since you have to actively erase
                      // it at the start of each loop
                      // And i BELIEVE that the "computational cost" of allocanting memory for a new variable
                      // is kinda the same as erasing a variable (i.e. overwriting it with "zeros")
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
