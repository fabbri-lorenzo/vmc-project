// Type definitions and constants

#ifndef VMCPROJECT_TYPES_HPP
#define VMCPROJECT_TYPES_HPP

// TODO: Use somewhere the Jackson-Freebeerg kinetic energy

#include <array>
#include <cassert>
#include <iterator>
#include <random>
#include <type_traits>

namespace vmcp {

// 'Structure' types
//
// Floating point type, can be adjusted to improve precision or compilation time
using FPType = long double;
// Unsigned integer type
// Only to be used for arrays, sizes etc.
// When it necessary to count something, use signed integers
using UIntType = long unsigned int;
// Signed integer type
using IntType = int;
// Tyoe for iterator’s difference type
using DiffType = std::vector<FPType>::difference_type;
// Used instead of default_random_engine since that one is implementation-defined
using RandomGenerator = std::mt19937;

// 'Program lexicon' types
//
// Dimension of the problem, usually it is 1, 2, 3
using Dimension = UIntType;
// Number of particles
using ParticNum = UIntType;
// Number of variational parameters
using VarParNum = UIntType;

// 'Wrapper' structs (so that it is impossible to add two objects of different types)
//
// Position of the particles (in D dimensions)
struct Coordinate {
    FPType val;
};
template <Dimension D>
using Position = std::array<Coordinate, D>;
template <Dimension D, ParticNum N>
using Positions = std::array<Position<D>, N>;
// Variational parameters
struct VarParam {
    FPType val;
};
template <VarParNum V>
using VarParams = std::array<VarParam, V>;
// Mass
struct Mass {
    FPType val;
};
// Variational Monte Carlo algorithm result
struct Energy {
    FPType val;
};
struct EnVariance {
    FPType val;
};
struct VMCResult {
    Energy energy;
    EnVariance variance;
};
struct BlockingResult {
    std::vector<IntType> sizes;
    std::vector<FPType> means;
    std::vector<FPType> stdDevs;
};
struct ConfInterval {
    FPType min;
    FPType max;
};
struct BootstrapResult {
    FPType meanOfMeans;
    FPType stdDevOfMeans;
    ConfInterval confInterval;
};
// TODO: Rename if the function is renamed
// WrappedVMVEnergies_ result
template <Dimension D, ParticNum N>
struct EnAndPos {
    Energy energy;
    Positions<D, N> positions;
};
// One-dimensional region of integration
struct Bound {
    FPType lower;
    FPType upper;
    Bound(FPType lower_, FPType upper_) : lower{lower_}, upper{upper_} { assert(upper >= lower); }
    FPType Length() const { return upper - lower; }
};
template <Dimension D>
using Bounds = std::array<Bound, D>;

// To (only) be used in a static assertion
template <Dimension D, ParticNum N, VarParNum V, class Function>
constexpr bool IsWavefunction() {
    return std::is_invocable_r_v<FPType, Function, Positions<D, N> const &, VarParams<V>>;
}
template <Dimension D, ParticNum N, VarParNum V, class Function>
constexpr bool IsWavefunctionDerivative() {
    return IsWavefunction<D, N, V, Function>();
}
template <Dimension D, ParticNum N, class Function>
constexpr bool IsPotential() {
    return std::is_invocable_r_v<FPType, Function, Positions<D, N> const &>;
}

} // namespace vmcp

#endif
