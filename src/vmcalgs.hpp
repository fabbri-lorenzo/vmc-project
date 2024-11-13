// Variational Monte Carlo algorithm definitions

#ifndef VMCPROJECT_VMCALGS_HPP
#define VMCPROJECT_VMCALGS_HPP

#include "types.hpp"

// TODO: Find a way to avoid having to always specify the template arguments (usually D, N), it's embarassing

namespace vmcp {

// Computes energies and positions by using the VMC algorithm with the analytical formula for the derivative
// and the Metropolis algorithm
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class SecondDerivative, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &, VarParams<V>, SecondDerivative const &,
                                                Mass, Potential const &, CoordBounds<D>, IntType,
                                                RandomGenerator &);

// Computes energy and variance by using the VMC algorithm with the analytical formula for the derivative and
// the Metropolis algorithm, after finding the best parameters
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class SecondDerivative, class Potential>
VMCResult VMCEnergy(Wavefunction const &, ParamBounds<V>, SecondDerivative const &, Mass, Potential const &,
                    CoordBounds<D>, IntType, RandomGenerator &);

// Computes energies and positions by using the VMC algorithm with the analytical formula for the derivative
// and the importance sampling algorithm
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative,
          class SecondDerivative, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &, VarParams<V>,
                                                std::array<FirstDerivative, D> const &,
                                                SecondDerivative const &, Mass, Potential const &,
                                                CoordBounds<D>, IntType, RandomGenerator &);

// Computes energy and variance by using the VMC algorithm with the analytical formula for the derivative and
// the importance sampling algorithm, after finding the best parameters
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative,
          class SecondDerivative, class Potential>
VMCResult VMCEnergy(Wavefunction const &, ParamBounds<V>, std::array<FirstDerivative, D> const &,
                    SecondDerivative const &, Mass, Potential const &, CoordBounds<D>, IntType,
                    RandomGenerator &);

// Computes energies and positions by using the VMC algorithm with a numerical estimation of the derivative
// and either the Metropolis or the importance sampling algorithm
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &, VarParams<V>, bool, FPType, Mass,
                                                Potential const &, CoordBounds<D>, IntType,
                                                RandomGenerator &);

// Computes energy and variance using the VMC algorithm with with a numerical estimation of the derivative
// and the either the Metropolis or the importance sampling algorithm, after finding the best parameters
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
VMCResult VMCEnergy(Wavefunction const &, ParamBounds<V>, bool, FPType, Mass, Potential const &,
                    CoordBounds<D>, IntType, RandomGenerator &);

} // namespace vmcp

// Implementation of templates is in this file
// It is separated to improve readability
#include "vmcalgs.inl"

#endif
