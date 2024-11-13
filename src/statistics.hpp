// Statistical analysis algorithms definitions

#ifndef VMCPROJECT_STATISTICS_HPP
#define VMCPROJECT_STATISTICS_HPP

#include "types.hpp"

namespace vmcp {
BlockingResult BlockingAnalysis(std::vector<Energy> const &energies);

BootstrapResult BootstrapAnalysis(std::vector<Energy> const &energies, UIntType const numSamples,
                                  RandomGenerator &gen);
} // namespace vmcp

// Implementation of templates is in this file
// It is separated to improve readability
#include "statistics.inl"

#endif
