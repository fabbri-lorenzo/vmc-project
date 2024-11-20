// Statistical analysis algorithms definitions

#ifndef VMCPROJECT_STATISTICS_HPP
#define VMCPROJECT_STATISTICS_HPP

#include "types.hpp"

namespace vmcp {

// Divides dataset into multiple blocks with fixed block size, then evaluate means of each block and
// carries on a statistical analysis of this vectors (one for each block size) of means
BlockingResult BlockingAnalysis(std::vector<Energy> const &energies);

// Samples dataset with replacement multiple times, then evaluates mean of each sample and
// carries on a statistical analysis of this vector of means
BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen);

} // namespace vmcp

// Implementation of templates is in this file
// It is separated to improve readability
#include "statistics.inl"

#endif
