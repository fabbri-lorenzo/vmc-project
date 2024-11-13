//
//
// Contains the definition of the templates declared in vmcalgs.hpp
// This file is supposed to only be #included at the end of vmcalgs.hpp and in vmcalgs.cpp
// It is just a way to improve the readability of vmcalgs.hpp
//

// The functions that have a trailing underscore would, in a regular project, be (declared and) defined in a
// .cpp file, therefore making them unreachable for the user
// Since they are templated, they must be defined in a header

#ifndef VMCPROJECT_VMCALGS_INL
#define VMCPROJECT_VMCALGS_INL

#include "vmcalgs.hpp"

#include <algorithm>
#include <cassert>
#include <functional>

namespace vmcp {

// Constants
// Their name is before the underscore
// After the underscore is the function(s) that use them
constexpr FPType hbar = 1; // For simplicity
constexpr IntType points_peakSearch = 100;
// LF TODO: Rename this one
// Also, is there a better criterion to choose this than to simply fix it (which is a bad idea since it might
// turn out to be too large for some programs)?
// I propose the following criterion: ask for 'step' as an input too (like MetropolisUpdate_ does), and choose
// deltaT such that on average the random part of the jump has length 'step'
constexpr FPType deltaT = 0.005f;
// FP TODO: Rename this one
constexpr IntType boundSteps = 100;
constexpr IntType thermalizationMoves = 100;
// FP TODO: Rename this one, also find if it is really 10 * thermalization moves
constexpr IntType movesForgetICs = 10 * thermalizationMoves;
constexpr FPType minPsi_peakSearch = 1e-6f;
constexpr IntType maxLoops_gradDesc = 10000;
constexpr IntType stepDenom_gradDesc = 100;
constexpr FPType stoppingThreshold_gradDesc = 1e-6f;
constexpr FPType targetAcceptRate_VMCLocEnAndPoss = 0.5f;
constexpr UIntType numWalkers_gradDesc = 10;
// FP TODO: Rename
constexpr FPType gradDescentFraction = 0.1f;

VMCResult AvgAndVar_(std::vector<Energy> const &);

template <Dimension D, ParticNum N>
std::vector<Energy> LocalEnergies_(std::vector<LocEnAndPoss<D, N>> const &v) {
    std::vector<Energy> result;
    result.reserve(v.size());
    std::transform(v.begin(), v.end(), std::back_inserter(result),
                   [](LocEnAndPoss<D, N> const &lep) { return lep.energy; });
    return result;
}

// Computes the position of the particles when the n-th one is moved by delta along the d-th direction (with
// both n and d starting from 0)
template <Dimension D, ParticNum N>
Positions<D, N> MoveBy_(Positions<D, N> const &poss, Dimension d, ParticNum n, Coordinate delta) {
    assert(d < D);
    assert(n < N);
    Positions<D, N> result = poss;
    result[n][d] += delta;
    return result;
}

// Randomly chooses points in the integration domain and computes where the potential is largest
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
Positions<D, N> FindPeak_(Wavefunction const &psi, VarParams<V> params, Potential const &pot,
                          CoordBounds<D> bounds, IntType points, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsPotential<D, N, Potential>());
    assert(points > 0);

    Position<D> center;
    // TODO: Define some operators for the coordinates?
    std::transform(bounds.begin(), bounds.end(), center.begin(),
                   [](Bound<Coordinate> b) { return Coordinate{(b.upper.val + b.lower.val) / 2}; });
    Positions<D, N> result;
    std::fill(result.begin(), result.end(), center);
    std::uniform_real_distribution<FPType> unif(0, 1);
    for (IntType i = 0; i != points; ++i) {
        Positions<D, N> newPoss;
        // FP TODO: Is it really necessary to specify Position<D>?
        // In general, study class template argument deduction
        for (Position<D> &p : newPoss) {
            std::transform(bounds.begin(), bounds.end(), p.begin(), [&unif, &gen](Bound<Coordinate> b) {
                return Coordinate{b.lower.val + unif(gen) * (b.upper.val - b.lower.val)};
            });
        }
        // The requirement ... > minPsi avoids having psi(...) = nan in the future, which breaks the update
        // algorithms
        if ((pot(newPoss) > pot(result)) && (psi(newPoss, params) > minPsi_peakSearch)) {
            result = newPoss;
        }
    }
    return result;
}

// Updates the wavefunction with the Metropolis algorithms and outputs the number of succesful updates
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction>
IntType MetropolisUpdate_(Wavefunction const &psi, VarParams<V> params, Positions<D, N> &poss, FPType step,
                          RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());

    IntType succesfulUpdates = 0;
    for (Position<D> &p : poss) {
        Position const oldPos = p;
        FPType const oldPsi = psi(poss, params);
        std::uniform_real_distribution<FPType> unif(0, 1);
        std::transform(p.begin(), p.end(), p.begin(), [&gen, &unif, step](Coordinate c) {
            return c + Coordinate{(unif(gen) - FPType{0.5f}) * step};
        });
        if (unif(gen) < std::pow(psi(poss, params) / oldPsi, 2)) {
            ++succesfulUpdates;
        } else {
            p = oldPos;
        }
    }
    return succesfulUpdates;
}

// Computes the drift force by using its analytic expression
template <Dimension D, ParticNum N, VarParNum V, class FirstDerivative, class Wavefunction>
std::array<FPType, D> DriftForceAnalytic_(Wavefunction const &psi, Positions<D, N> poss, VarParams<V> params,
                                          std::array<FirstDerivative, D> const &grad) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsWavefunctionDerivative<D, N, V, FirstDerivative>());

    std::array<FPType, D> gradVal;
    std::transform(grad.begin(), grad.end(), gradVal.begin(),
                   [&poss, params](FirstDerivative const &fd) { return fd(poss, params); });
    std::array<FPType, D> result;
    std::transform(gradVal.begin(), gradVal.end(), result.begin(),
                   [&psi, &poss, params](FPType f) { return 2 * f / psi(poss, params); });
    return result;
}

// Computes the drift force by numerically estimating the derivative of the wavefunction
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction>
std::array<FPType, D> DriftForceNumeric_(Wavefunction const &psi, VarParams<V> params, FPType derivativeStep,
                                         Positions<D, N> poss) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());

    std::array<FPType, D> driftForce;
    for (ParticNum n = 0u; n != N; ++n) {
        for (Dimension d = 0u; d != D; ++d) {
            driftForce[d] = 2 *
                            (psi(MoveBy_<D, N>(poss, d, n, Coordinate{derivativeStep}), params) -
                             psi(MoveBy_<D, N>(poss, d, n, Coordinate{-derivativeStep}), params)) /
                            (derivativeStep * psi(poss, params));
        }
    }
}

// LF TODO: Tell that this functions applies formula at page 22 Jensen
// LF TODO: DriftForceNumeric_ is not used here!
//  Updates the wavefunction with the importance sampling algorithm and outputs the number of succesful
//  updates
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative>
IntType ImportanceSamplingUpdate_(Wavefunction const &psi, VarParams<V> params,
                                  std::array<FirstDerivative, D> const &grad, Mass mass,
                                  Positions<D, N> &poss, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsWavefunctionDerivative<D, N, V, FirstDerivative>());

    FPType const diffusionConst = hbar * hbar / (2 * mass.val);

    IntType successfulUpdates = 0;
    for (Position<D> &p : poss) {
        Position const oldPos = p;
        FPType const oldPsi = psi(poss, params);
        std::array<FPType, D> const oldDriftForce = DriftForceAnalytic_<D, N, V>(psi, poss, params, grad);
        std::normal_distribution<FPType> normal(0, diffusionConst * deltaT);
        for (Dimension d = 0u; d != D; ++d) {
            p[d].val = oldPos[d].val + diffusionConst * deltaT * (oldDriftForce[d] + normal(gen));
        }
        FPType const newPsi = psi(poss, params);
        std::array<FPType, D> const newDriftForce = DriftForceAnalytic_<D, N, V>(psi, poss, params, grad);
        FPType forwardExponent = 0;
        for (Dimension d = 0u; d != D; ++d) {
            forwardExponent -=
                std::pow(p[d].val - oldPos[d].val - diffusionConst * deltaT * oldDriftForce[d], 2) /
                (4 * diffusionConst * deltaT);
        }
        FPType const forwardProb = std::exp(forwardExponent);
        FPType backwardExponent = 0;
        for (Dimension d = 0u; d != D; ++d) {
            backwardExponent -=
                std::pow(oldPos[d].val - p[d].val - diffusionConst * deltaT * newDriftForce[d], 2) /
                (4 * diffusionConst * deltaT);
        }
        FPType const backwardProb = std::exp(backwardExponent);
        FPType const acceptanceRatio = (newPsi * newPsi * backwardProb) / (oldPsi * oldPsi * forwardProb);
        std::uniform_real_distribution<FPType> unif(0, 1);
        if (unif(gen) < acceptanceRatio) {
            ++successfulUpdates;
        } else {
            p = oldPos;
        }
    }
    return successfulUpdates;
}

// Computes the local energy with the analytic formula for the derivative of the wavefunction
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class SecondDerivative, class Potential>
Energy LocalEnergyAnalytic_(Wavefunction const &psi, VarParams<V> params, SecondDerivative const &secondDer,
                            Mass mass, Potential const &pot, Positions<D, N> poss) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsWavefunctionDerivative<D, N, V, SecondDerivative>());
    static_assert(IsPotential<D, N, Potential>());

    return Energy{-(hbar * hbar / (2 * mass.val)) * (secondDer(poss, params) / psi(poss, params)) +
                  pot(poss)};
}

// Computes the local energy by numerically estimating the derivative of the wavefunction
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
Energy LocalEnergyNumeric_(Wavefunction const &psi, VarParams<V> params, FPType derivativeStep, Mass mass,
                           Potential const &pot, Positions<D, N> poss) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsPotential<D, N, Potential>());

    Energy result{pot(poss)};
    for (ParticNum n = 0u; n != N; ++n) {
        for (Dimension d = 0u; d != D; ++d) {
            result.val +=
                -std::pow(hbar, 2) / (2 * mass.val) *
                (psi(MoveBy_<D, N>(poss, d, n, Coordinate{derivativeStep}), params) - 2 * psi(poss, params) +
                 psi(MoveBy_<D, N>(poss, d, n, Coordinate{-derivativeStep}), params)) /
                std::pow(derivativeStep, 2);
        }
    }
    return result;
}

// Computes the energies that will be averaged to obtain the estimate of the GS energy of the VMC algorithm
// by either using the analytical formula for the derivative or estimating it numerically
// Use the wrappers to select which method should be used
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative,
          class SecondDerivative, class Potential>
std::vector<LocEnAndPoss<D, N>>
VMCLocEnAndPoss_(Wavefunction const &psi, VarParams<V> params, bool useAnalytical, bool useImpSamp,
                 std::array<FirstDerivative, D> const &grad, SecondDerivative const &secondDer,
                 FPType derivativeStep, Mass mass, Potential const &pot, CoordBounds<D> bounds,
                 IntType numberEnergies, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsWavefunctionDerivative<D, N, V, FirstDerivative>());
    static_assert(IsWavefunctionDerivative<D, N, V, SecondDerivative>());
    static_assert(IsPotential<D, N, Potential>());
    assert(numberEnergies > 0);

    std::vector<LocEnAndPoss<D, N>> result;

    // Find a good starting point, in the sense that it's easy to move away from
    Positions<D, N> const peak = FindPeak_<D, N>(psi, params, pot, bounds, points_peakSearch, gen);

    // Choose the step
    Bound const smallestBound =
        *(std::min_element(bounds.begin(), bounds.end(), [](Bound<Coordinate> b1, Bound<Coordinate> b2) {
            return b1.Length().val < b2.Length().val;
        }));
    FPType step = smallestBound.Length().val / boundSteps;

    // Save the energies to be averaged
    Positions<D, N> poss = peak;
    std::function<Energy()> localEnergy;
    if (useAnalytical) {
        localEnergy = std::function<Energy()>{
            [&]() { return LocalEnergyAnalytic_<D, N>(psi, params, secondDer, mass, pot, poss); }};
    } else {
        localEnergy = std::function<Energy()>{
            [&]() { return LocalEnergyNumeric_<D, N>(psi, params, derivativeStep, mass, pot, poss); }};
    }
    std::function<IntType()> update;
    if (useImpSamp) {
        update = std::function<IntType()>{
            [&]() { return ImportanceSamplingUpdate_<D, N>(psi, params, grad, mass, poss, gen); }};
    } else {
        update =
            std::function<IntType()>{[&]() { return MetropolisUpdate_<D, N>(psi, params, poss, step, gen); }};
    }

    // Move away from the peak, in order to forget the dependence on the initial conditions
    for (IntType i = 0; i != movesForgetICs; ++i) {
        update();
    }
    for (IntType i = 0; i != numberEnergies; ++i) {
        IntType succesfulUpdates = 0;
        for (IntType j = 0; j != thermalizationMoves; ++j) {
            succesfulUpdates += update();
        }
        result.emplace_back(localEnergy(), poss);

        // Adjust the step size
        // Call car = current acc. rate, tar = target acc. rate
        // Add (car - tar)/tar to step, since it increases step if too many moves were accepted and decreases
        // it if too few were accepted
        FPType currentAcceptRate = succesfulUpdates * FPType{1} / (thermalizationMoves * N);
        step *= currentAcceptRate / targetAcceptRate_VMCLocEnAndPoss;
    }

    return result;
}

// Computes the energies that would be calculated by the VMC algorithm after the parameters are moved in one
// cardinal direction, by using the reweighting method
// Used to compute the gradient of the VMC energy in parameter space
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction>
std::array<Energy, D> ReweightedEnergies_(Wavefunction const &psi, VarParams<V> oldParams,
                                          std::vector<LocEnAndPoss<N, D>> oldLEPs, FPType step) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());

    std::array<Energy, V> result;
    for (UIntType v = 0u; v != V; ++v) {
        VarParams<V> newParams = oldParams;
        newParams[v].val += step;
        std::vector<Energy> reweightedLocalEnergies;
        std::transform(oldLEPs.begin(), oldLEPs.end(), std::back_inserter(reweightedLocalEnergies),
                       [&psi, newParams, oldParams](LocEnAndPoss<D, N> const &lep) {
                           return Energy{
                               std::pow(psi(lep.positions, newParams) / psi(lep.positions, oldParams), 2) *
                               lep.energy.val};
                       });

        FPType const numerator =
            std::accumulate(reweightedLocalEnergies.begin(), reweightedLocalEnergies.end(), FPType{0},
                            [](FPType f, Energy e) { return f + e.val; });
        FPType const denominator = std::accumulate(
            oldLEPs.begin(), oldLEPs.end(), FPType{0},
            [&psi, newParams, oldParams](FPType f, LocEnAndPoss<D, N> const &lep) {
                return f + std::pow(psi(lep.positions, newParams) / psi(lep.positions, oldParams), 2);
            });

        result[v] = Energy{numerator / denominator};
    }
    return result;
}

// Calculates the parameters of the wavefunction that minimize the VMC energy, by doing gradient descent
// starting from the given parameters
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class LocEnAndPossCalculator>
VMCResult VMCRBestParams_(VarParams<V> initialParams, Wavefunction const &psi,
                          LocEnAndPossCalculator const &lepsCalc) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(
        std::is_invocable_r_v<std::vector<LocEnAndPoss<D, N>>, LocEnAndPossCalculator, VarParams<V>>);
    static_assert(V != UIntType{0});

    // Use a modified gradient descent algorithm with termination condition: stop if the new parameters
    // are too close to the old ones
    VMCResult result;
    VarParams<V> currentParams = initialParams;
    FPType const initialParamsNorm = std::sqrt(
        std::inner_product(initialParams.begin(), initialParams.end(), initialParams.begin(), FPType{0},
                           std::plus<>(), [](VarParam v1, VarParam v2) { return v1.val * v2.val; }));
    FPType gradientStep = initialParamsNorm / stepDenom_gradDesc;

    for (IntType i = 0; i != maxLoops_gradDesc; ++i) {
        std::vector<LocEnAndPoss<D, N>> const currentEnAndPoss = lepsCalc(currentParams);
        VMCResult const currentVMCR = AvgAndVar_(LocalEnergies_(currentEnAndPoss));
        result = currentVMCR;

        // Compute the gradient by using reweighting
        std::array<Energy, V> energiesIncreasedParam =
            ReweightedEnergies_<D, N, V>(psi, currentParams, currentEnAndPoss, gradientStep);
        std::array<Energy, V> energiesDecreasedParam =
            ReweightedEnergies_<D, N, V>(psi, currentParams, currentEnAndPoss, -gradientStep);
        std::array<FPType, V> gradient;
        for (VarParNum v = 0u; v != V; ++v) {
            gradient[v] =
                (energiesIncreasedParam[v].val - energiesDecreasedParam[v].val) / (2 * gradientStep);
        }

        // Set as next step used to compute the gradient the current gradient norm, which is also the size of
        // the step if that step is accepted
        FPType const oldParamsNorm =
            std::sqrt(std::accumulate(currentParams.begin(), currentParams.end(), FPType{0},
                                      [](FPType f, VarParam v) { return f + v.val * v.val; }));
        gradientStep =
            std::sqrt(std::inner_product(gradient.begin(), gradient.end(), gradient.begin(), FPType{0}));
        if ((gradientStep / oldParamsNorm) < stoppingThreshold_gradDesc) {
            break;
        }

        // In order to move, try smaller and smaller steps until the new parameters produce a VMC energy not
        // much larger than the current one
        FPType gradMultiplier = 2;
        VarParams<V> newParams;
        VMCResult newVMCR;
        do {
            gradMultiplier /= 2;
            for (VarParNum v = 0u; v != V; ++v) {
                newParams[v].val = currentParams[v].val - gradMultiplier * gradient[v];
                newVMCR = AvgAndVar_(LocalEnergies_(lepsCalc(newParams)));
            }

        } while (newVMCR.energy.val >
                 (currentVMCR.energy.val + gradDescentFraction * std::abs(currentVMCR.energy.val)));
        currentParams = newParams;
    }

    return result;
}

// Calculates the parameters of the wavefunction that minimize the VMC energy, by doing gradient descent from
// many different initial starting points and picking the one that arrived to the lowest point
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class LocEnAndPossCalculator>
VMCResult VMCRBestParams_(ParamBounds<V> bounds, Wavefunction const &psi,
                          LocEnAndPossCalculator const &lepsCalc, UIntType numWalkers, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(
        std::is_invocable_r_v<std::vector<LocEnAndPoss<D, N>>, LocEnAndPossCalculator, VarParams<V>>);

    if constexpr (V == VarParNum{0u}) {
        VarParams<0> const fakeParams{};
        return AvgAndVar_(LocalEnergies_(lepsCalc(fakeParams)));
    } else {
        std::uniform_real_distribution<FPType> unif(0, 1);
        std::vector<VMCResult> vmcResults;
        vmcResults.reserve(numWalkers);
        std::generate_n(std::back_inserter(vmcResults), numWalkers_gradDesc, [&]() {
            VarParams<V> initialParams;
            for (VarParNum v = 0u; v != V; ++v) {
                initialParams[v].val = bounds[v].lower.val + bounds[v].Length().val * unif(gen);
            }
            return VMCRBestParams_<D, N, V>(initialParams, psi, lepsCalc);
        });
        return *std::min_element(vmcResults.begin(), vmcResults.end(), [](VMCResult vmcr1, VMCResult vmcr2) {
            return vmcr1.energy.val < vmcr2.energy.val;
        });
    }
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class SecondDerivative, class Potential>
std::vector<LocEnAndPoss<D, N>>
VMCLocEnAndPoss(Wavefunction const &psi, VarParams<V> params, SecondDerivative const &secondDer, Mass mass,
                Potential const &pot, CoordBounds<D> bounds, IntType numberEnergies, RandomGenerator &gen) {
    struct FakeDeriv {
        FPType operator()(Positions<D, N> const &, VarParams<V>) const {
            assert(false);
            return FPType{0};
        }
    };
    std::array<FakeDeriv, D> fakeGrad;
    std::fill(fakeGrad.begin(), fakeGrad.end(), FakeDeriv{});
    FPType const fakeStep = std::numeric_limits<FPType>::quiet_NaN();
    return VMCLocEnAndPoss_<D, N, V>(psi, params, true, false, fakeGrad, secondDer, fakeStep, mass, pot,
                                     bounds, numberEnergies, gen);
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class SecondDerivative, class Potential>
VMCResult VMCEnergy(Wavefunction const &psi, ParamBounds<V> parBounds, SecondDerivative const &secondDer,
                    Mass mass, Potential const &pot, CoordBounds<D> coorBounds, IntType numberEnergies,
                    RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(psi, vps, secondDer, mass, pot, coorBounds, numberEnergies, gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, psi, enPossCalculator, numWalkers_gradDesc, gen);
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative,
          class SecondDerivative, class Potential>
std::vector<LocEnAndPoss<D, N>>
VMCLocEnAndPoss(Wavefunction const &psi, VarParams<V> params, std::array<FirstDerivative, D> const &grad,
                SecondDerivative const &secondDer, Mass mass, Potential const &pot, CoordBounds<D> bounds,
                IntType numberEnergies, RandomGenerator &gen) {
    FPType const fakeStep = std::numeric_limits<FPType>::quiet_NaN();
    return VMCLocEnAndPoss_<D, N, V>(psi, params, true, false, grad, secondDer, fakeStep, mass, pot, bounds,
                                     numberEnergies, gen);
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative,
          class SecondDerivative, class Potential>
VMCResult VMCEnergy(Wavefunction const &psi, ParamBounds<V> parBounds,
                    std::array<FirstDerivative, D> const &grad, SecondDerivative const &secondDer, Mass mass,
                    Potential const &pot, CoordBounds<D> coorBounds, IntType numberEnergies,
                    RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(psi, vps, grad, secondDer, mass, pot, coorBounds, numberEnergies,
                                        gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, psi, enPossCalculator, numWalkers_gradDesc, gen);
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &psi, VarParams<V> params, bool useImpSamp,
                                                FPType derivativeStep, Mass mass, Potential const &pot,
                                                CoordBounds<D> bounds, IntType numberEnergies,
                                                RandomGenerator &gen) {
    struct FakeDeriv {
        FPType operator()(Positions<D, N> const &, VarParams<V>) const {
            assert(false);
            return FPType{0};
        }
    };
    std::array<FakeDeriv, D> fakeGrad;
    std::fill(fakeGrad.begin(), fakeGrad.end(), FakeDeriv{});
    return VMCLocEnAndPoss_<D, N, V>(psi, params, false, useImpSamp, fakeGrad, FakeDeriv{}, derivativeStep,
                                     mass, pot, bounds, numberEnergies, gen);
}

template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
VMCResult VMCEnergy(Wavefunction const &psi, ParamBounds<V> parBounds, bool useImpSamp, FPType derivativeStep,
                    Mass mass, Potential const &pot, CoordBounds<D> coorBounds, IntType numberEnergies,
                    RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(psi, vps, useImpSamp, derivativeStep, mass, pot, coorBounds,
                                        numberEnergies, gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, psi, enPossCalculator, numWalkers_gradDesc, gen);
}

} // namespace vmcp

#endif
