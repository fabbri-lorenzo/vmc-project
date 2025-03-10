//!
//! @file vmcalgs.inl
//! @brief Definition of the user and core functions
//! @authors Lorenzo Fabbri, Francesco Orso Pancaldi
//!
//! The user functions are declared in the related .hpp
//! To improve readability, the definitions of the helper functions is in the other .inl
//! @see vmcalgs.hpp
//! @see vmchelpers.inl
//!

#ifndef VMCPROJECT_VMCALGS_INL
#define VMCPROJECT_VMCALGS_INL

#include "statistics.hpp"
#include "vmcalgs.hpp"
#include "vmchelpers.inl"

namespace vmcp {

//! @defgroup core-functions Core functions
//! @brief The most important functions in the code
//!
//! The ones that actually do the work.
//! @{

//! @brief Computes the energies that will be averaged to obtain the estimate of the energy
//! @param wavef The wavefunction
//! @param params The variational parameters
//! @param useAnalytical Whether the local energy must be computed by using the analytical expression of the
//! laplacians
//! @param useImpSamp Whether to use importance sampling as the update algorithm (the alternative is
//! Metropolis)
//! @param grads The gradients of the particles (unused if 'useAnalytical == false' or 'useImpSamp == false')
//! @param lapls The laplacians of the particles (unused if 'useAnalytical == false')
//! @param derivativeStep The step used is the numerical estimation of the derivative (unused if
//! 'useAnalytical == true')
//! @param masses The masses of the particles
//! @param pot The potential
//! @param bounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Is the most important function of the library, since is the one which actually does the work.
//! Starts from a point where the potential is sufficiently large, to quickly forget about the initial
//! conditions. Does some updates to move away from the starting point, then starts computing the local
//! energies. In between two evaluations of the local energy, some updates are done to avoid correlations.
//! Adjusts the step size on the fly to best match the target acceptance rate. Depending on 'useAnalytical'
//! and 'useImpSamp', some parameters are unused. To avoid having the user supply some parameters he does not
//! care about, wrappers that only ask for the necessary ones are provided.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative, class Laplacian,
          class Potential>
std::vector<LocEnAndPoss<D, N>>
VMCLocEnAndPoss_(Wavefunction const &wavef, Positions<D, N> poss, VarParams<V> params, bool useAnalytical,
                 bool useImpSamp, Gradients<D, N, FirstDerivative> const &grads,
                 Laplacians<N, Laplacian> const &lapls, FPType derivativeStep, Masses<N> masses,
                 Potential const &pot, CoordBounds<D> bounds, IntType numEnergies, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(IsWavefunctionDerivative<D, N, V, FirstDerivative>());
    static_assert(IsWavefunctionDerivative<D, N, V, Laplacian>());
    static_assert(IsPotential<D, N, Potential>());
    assert(numEnergies > 0);

    // Choose the initial step
    Bound const smallestBound =
        *(std::min_element(bounds.begin(), bounds.end(), [](Bound<Coordinate> b1, Bound<Coordinate> b2) {
            return b1.Length().val < b2.Length().val;
        }));
    FPType step = smallestBound.Length().val / stepDenom_vmcLEPs;

    std::function<Energy()> localEnergy;
    if (useAnalytical) {
        localEnergy = std::function<Energy()>{
            [&]() { return LocalEnergyAnalytic_<D, N>(wavef, params, lapls, masses, pot, poss); }};
    } else {
        localEnergy = std::function<Energy()>{
            [&]() { return LocalEnergyNumeric_<D, N>(wavef, params, derivativeStep, masses, pot, poss); }};
    }
    std::function<IntType()> update;
    if (useImpSamp) {
        update = std::function<IntType()>{[&]() {
            return ImportanceSamplingUpdate_<D, N>(wavef, params, useAnalytical, derivativeStep, grads,
                                                   masses, poss, gen);
        }};
    } else {
        update = std::function<IntType()>{
            [&]() { return MetropolisUpdate_<D, N>(wavef, params, poss, step, gen); }};
    }

    std::vector<LocEnAndPoss<D, N>> result;
    result.reserve(static_cast<long unsigned int>(numEnergies));
    // Move away from the starting ponit, in order to forget the dependence on the initial conditions
    for (IntType i = 0; i != movesForgetICs_vmcLEPs; ++i) {
        update();
    }
    for (IntType i = 0; i != numEnergies; ++i) {
        IntType succesfulUpdates = 0;
        for (IntType j = 0; j != autocorrelationMoves_vmcLEPs; ++j) {
            succesfulUpdates += update();
        }
        result.emplace_back(localEnergy(), poss);

        // Adjust the step size
        // Call car = current acc. rate, tar = target acc. rate
        // Add (car - tar)/tar to step, since it increases step if too many moves were accepted and decreases
        // it if too few were accepted
        FPType currentAcceptRate = succesfulUpdates * FPType{1} / (autocorrelationMoves_vmcLEPs * N);
        step *= (currentAcceptRate > targetAcceptRate_vmcLEPs ? FPType{11} / 10 : FPType{9} / 10);
    }

    return result;
}

//! @brief Computes the energy with error, with the parameters that minimize the former
//! @param initialParams The initial variational parameters
//! @param wavef The wavefunction
//! @param lepsCalc A function that takes as input the variational parameters and returns the local energies
//! and the positions of the particles when each one was computed
//! @return The energy with error
//!
//! Does gradient descent starting from the given parameters.
//! Stops when the proposed step is too small compared to the current parameters.
//! Computes the gradient by using reweighting.
//! After having computed the gradient, if the proposed step would increase the energy too much, proposes a
//! new step of half the length, and repeats.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class LocEnAndPossCalculator>
VMCResult<V> VMCRBestParams_(VarParams<V> initialParams, ParamBounds<V> bounds, Wavefunction const &wavef,
                             LocEnAndPossCalculator const &lepsCalc, StatFuncType function,
                             IntType const &boostrapSamples, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(
        std::is_invocable_r_v<std::vector<LocEnAndPoss<D, N>>, LocEnAndPossCalculator, VarParams<V>>);
    static_assert(V != 0);
    assert(!std::isnan(initialParams[0].val));

    // Use a gradient descent algorithm with termination condition: stop if the next step is too small
    // compared to the current parameters
    VMCResult<V> result;
    Energy currentEn;
    VarParams<V> currentParams = initialParams;
    FPType const initialParamsNorm = std::sqrt(
        std::inner_product(initialParams.begin(), initialParams.end(), initialParams.begin(), FPType{0},
                           std::plus<>(), [](VarParam v1, VarParam v2) { return v1.val * v2.val; }));
    FPType gradStep = initialParamsNorm / stepDenom_gradDesc;
    std::array<FPType, V> oldMomentum;
    std::fill_n(oldMomentum.begin(), V, FPType{0});

    for (IntType i = 0; i != maxLoops_gradDesc; ++i) {
        // The gradient descent should end in a reasonable time
        assert((i + 1) != maxLoops_gradDesc);
        // Better to be turned off if numWalkers_gradDesc != 1
        std::cout << "Variational Parameter: " << currentParams[0].val << "\n";

        // Update the energy
        std::vector<LocEnAndPoss<D, N>> const currentLEPs = lepsCalc(currentParams);
        currentEn = Mean(currentLEPs);

        // Compute the gradient by using reweighting and update the momentum
        std::array<Energy, V> energiesIncreasedParam =
            ReweightedEnergies_<D, N, V>(wavef, currentParams, currentLEPs, gradStep);
        std::array<Energy, V> energiesDecreasedParam =
            ReweightedEnergies_<D, N, V>(wavef, currentParams, currentLEPs, -gradStep);
        std::array<FPType, V> currentMomentum;
        std::generate_n(currentMomentum.begin(), V,
                        [v = VarParNum{0}, &energiesIncreasedParam, &energiesDecreasedParam, &oldMomentum,
                         gradStep]() mutable {
                            FPType const result_ =
                                -FPType{3} / 4 *
                                    (energiesIncreasedParam[v].val - energiesDecreasedParam[v].val) /
                                    (2 * gradStep) +
                                FPType{1} / 4 * oldMomentum[v];
                            assert(!std::isnan(result_));
                            ++v;
                            return result_;
                        });

        // Set as next step used to compute the gradient the current gradient norm, which is also the size of
        // the step if that step is accepted
        FPType const currentParamsNorm = std::sqrt(
            std::inner_product(currentParams.begin(), currentParams.end(), currentParams.begin(), FPType{0},
                               std::plus<>(), [](VarParam v1, VarParam v2) { return v1.val * v2.val; }));
        gradStep = std::sqrt(std::inner_product(currentMomentum.begin(), currentMomentum.end(),
                                                currentMomentum.begin(), FPType{0}));

        // Check the termination condition
        if (gradStep / currentParamsNorm < stoppingThreshold_gradDesc) {
            result.energy = currentEn;
            result.stdDev = ErrorOnAvg(currentLEPs, function, boostrapSamples, gen);
            result.bestParams = currentParams;
            break;
        } else {
            for (VarParNum v = 0u; v != V; ++v) {
                FPType multiplier = 0.02f;
                while ((currentParams[v].val + multiplier * currentMomentum[v] > bounds[v].upper.val) ||
                       (currentParams[v].val + multiplier * currentMomentum[v] < bounds[v].lower.val)) {
                    multiplier /= 2;
                }
                currentParams[v].val += multiplier * currentMomentum[v];
            }
            /*std::transform(currentParams.begin(), currentParams.end(), currentMomentum.begin(),
                           currentParams.begin(), [](VarParam cp, FPType cm) { return cp + VarParam{cm}; });*/
            oldMomentum = currentMomentum;
        }
    }

    return result;
}

//! @brief Computes the energy with error, with the parameters that minimize the former
//! @param bounds The interval in which the best parameters should be found
//! @param wavef The wavefunction
//! @param lepsCalc A function that takes as input the variational parameters and returns the local energies
//! and the positions of the particles when each one was computed
//! @param numWalkers The number of independent gradient descents carried out
//! @param gen The random generator
//! @return The energy with error
//!
//! Carries out 'numWalkers' gradient descents in parallel, and at the end chooses the lowest energy obtained.
//! The starting parameters of the walkers are chosen randomly inside 'bounds'
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class LocEnAndPossCalculator>
VMCResult<V> VMCRBestParams_(ParamBounds<V> bounds, Wavefunction const &wavef,
                             LocEnAndPossCalculator const &lepsCalc, IntType numWalkers,
                             StatFuncType function, IntType const &boostrapSamples, RandomGenerator &gen) {
    static_assert(IsWavefunction<D, N, V, Wavefunction>());
    static_assert(
        std::is_invocable_r_v<std::vector<LocEnAndPoss<D, N>>, LocEnAndPossCalculator, VarParams<V>>);
    assert(numWalkers > IntType{0});

    if constexpr (V == VarParNum{0}) {
        VarParams<0u> const fakeParams{};
        std::vector<LocEnAndPoss<D, N>> const vmcLEPs = lepsCalc(fakeParams);
        return VMCResult<0>{Mean(vmcLEPs), ErrorOnAvg(vmcLEPs, function, boostrapSamples, gen),
                            VarParams<0>{}};
    } else {
        std::mutex m;
        std::uniform_real_distribution<FPType> unif(0, 1);
        std::vector<VMCResult<V>> vmcResults(static_cast<long unsigned int>(numWalkers));
        std::generate_n(std::execution::par_unseq, vmcResults.begin(), numWalkers_gradDesc, [&]() {
            unsigned long int seed_;
            {
                std::lock_guard<std::mutex> l(m);
                seed_ = gen();
            }
            RandomGenerator localGen{seed_};
            VarParams<V> initialParams;
            for (VarParNum v = 0u; v != V; ++v) {
                initialParams[v] = bounds[v].lower + bounds[v].Length() * unif(localGen);
            }
            return VMCRBestParams_<D, N, V>(initialParams, bounds, wavef, lepsCalc, function, boostrapSamples,
                                            localGen);
        });

        return *std::min_element(
            vmcResults.begin(), vmcResults.end(),
            [](VMCResult<V> vmcr1, VMCResult<V> vmcr2) { return vmcr1.energy < vmcr2.energy; });
    }
}

//! @}

//! @defgroup user-functions User functions
//! @brief The functions that are meant to be called by the user
//!
//! Are wrappers for the core functions.
//! @{

//! @brief Computes the energies that will be averaged by using the analytical formula for the derivative and
//! the Metropolis algorithm
//! @param wavef The wavefunction
//! @param params The variational parameters
//! @param lapls The laplacians of the particles
//! @param masses The masses of the particles
//! @param pot The potential
//! @param bounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for the true 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Laplacian, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &wavef, Positions<D, N> poss,
                                                VarParams<V> params, Laplacians<N, Laplacian> const &lapls,
                                                Masses<N> masses, Potential const &pot, CoordBounds<D> bounds,
                                                IntType numEnergies, RandomGenerator &gen) {
    struct FakeDeriv {
        FPType operator()(Positions<D, N> const &, VarParams<V>) const {
            assert(false);
            return 0;
        }
    };
    Gradients<D, N, FakeDeriv> fakeGrads;
    FPType const fakeStep = std::numeric_limits<FPType>::quiet_NaN();
    return VMCLocEnAndPoss_<D, N, V>(wavef, poss, params, true, false, fakeGrads, lapls, fakeStep, masses,
                                     pot, bounds, numEnergies, gen);
}

//! @brief Computes the energy with error, by using the analytical formula for the derivative and the
//! Metropolis algorithm, after finding the best parameter
//! @param wavef The wavefunction
//! @param parBounds The interval in which the best parameters should be found
//! @param lapls The laplacians of the particles
//! @param masses The masses of the particles
//! @param pot The potential
//! @param coorBounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Laplacian, class Potential>
VMCResult<V> VMCEnergy(Wavefunction const &wavef, Positions<D, N> poss, ParamBounds<V> parBounds,
                       Laplacians<N, Laplacian> const &lapls, Masses<N> masses, Potential const &pot,
                       CoordBounds<D> coorBounds, IntType numEnergies, StatFuncType function,
                       IntType const &boostrapSamples, RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(wavef, poss, vps, lapls, masses, pot, coorBounds, numEnergies, gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, wavef, enPossCalculator, numWalkers_gradDesc, function,
                                    boostrapSamples, gen);
}

//! @brief Computes the energies that will be averaged by using the analytical formula for the derivative and
//! the importance sampling algorithm
//! @param wavef The wavefunction
//! @param params The variational parameters
//! @param grads The gradients of the particles
//! @param lapls The laplacians of the particles
//! @param masses The masses of the particles
//! @param pot The potential
//! @param bounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for the true 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative, class Laplacian,
          class Potential>
std::vector<LocEnAndPoss<D, N>>
VMCLocEnAndPoss(Wavefunction const &wavef, Positions<D, N> poss, VarParams<V> params,
                Gradients<D, N, FirstDerivative> const &grads, Laplacians<N, Laplacian> const &lapls,
                Masses<N> masses, Potential const &pot, CoordBounds<D> bounds, IntType numEnergies,
                RandomGenerator &gen) {
    FPType const fakeStep = std::numeric_limits<FPType>::quiet_NaN();
    return VMCLocEnAndPoss_<D, N, V>(wavef, poss, params, true, true, grads, lapls, fakeStep, masses, pot,
                                     bounds, numEnergies, gen);
}

//! @brief Computes the energy with error, by using the analytical formula for the derivative and the
//! importance sampling algorithm, after finding the best parameter
//! @param wavef The wavefunction
//! @param parBounds The interval in which the best parameters should be found
//! @param grads The gradients of the particles
//! @param lapls The laplacians of the particles
//! @param masses The masses of the particles
//! @param pot The potential
//! @param coorBounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class FirstDerivative, class Laplacian,
          class Potential>
VMCResult<V> VMCEnergy(Wavefunction const &wavef, Positions<D, N> poss, ParamBounds<V> parBounds,
                       Gradients<D, N, FirstDerivative> const &grads, Laplacians<N, Laplacian> const &lapls,
                       Masses<N> masses, Potential const &pot, CoordBounds<D> coorBounds, IntType numEnergies,
                       StatFuncType function, IntType const &boostrapSamples, RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(wavef, poss, vps, grads, lapls, masses, pot, coorBounds, numEnergies,
                                        gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, wavef, enPossCalculator, numWalkers_gradDesc, function,
                                    boostrapSamples, gen);
}

//! @brief Computes the energies that will be averaged by numerically estimating the derivative and using
//! either the Metropolis or the importance sampling algorithm
//! @param wavef The wavefunction
//! @param params The variational parameters
//! @param useImpSamp Whether to use importance sampling as the update algorithm (the alternative is
//! Metropolis)
//! @param derivativeStep The step used is the numerical estimation of the derivative
//! @param masses The masses of the particles
//! @param pot The potential
//! @param bounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for the true 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
std::vector<LocEnAndPoss<D, N>> VMCLocEnAndPoss(Wavefunction const &wavef, Positions<D, N> poss,
                                                VarParams<V> params, bool useImpSamp, FPType derivativeStep,
                                                Masses<N> masses, Potential const &pot, CoordBounds<D> bounds,
                                                IntType numEnergies, RandomGenerator &gen) {
    struct FakeDeriv {
        FPType operator()(Positions<D, N> const &, VarParams<V>) const {
            assert(false);
            return FPType{0};
        }
    };
    Gradients<D, N, FakeDeriv> fakeGrads;
    std::array<FakeDeriv, N> fakeLapls;
    return VMCLocEnAndPoss_<D, N, V>(wavef, poss, params, false, useImpSamp, fakeGrads, fakeLapls,
                                     derivativeStep, masses, pot, bounds, numEnergies, gen);
}

//! @brief Computes the energy with error, by numerically estimating the derivative and using either the
//! Metropolis or the importance sampling algorithm, after finding the best parameter
//! @param wavef The wavefunction
//! @param parBounds The interval in which the best parameters should be found
//! @param useImpSamp Whether to use importance sampling as the update algorithm (the alternative is
//! Metropolis)
//! @param derivativeStep The step used in the numerical estimation of the derivative(s)
//! @param masses The masses of the particles
//! @param pot The potential
//! @param coorBounds The integration region
//! @param numEnergies The number of energies to compute
//! @param gen The random generator
//! @return The computed local energies, and the positions of the particles when each local energy was
//! computed
//!
//! Wrapper for 'VMCLocEnAndPoss'.
template <Dimension D, ParticNum N, VarParNum V, class Wavefunction, class Potential>
VMCResult<V> VMCEnergy(Wavefunction const &wavef, Positions<D, N> poss, ParamBounds<V> parBounds,
                       bool useImpSamp, FPType derivativeStep, Masses<N> masses, Potential const &pot,
                       CoordBounds<D> coorBounds, IntType numEnergies, StatFuncType function,
                       IntType const &boostrapSamples, RandomGenerator &gen) {
    auto const enPossCalculator{[&](VarParams<V> vps) {
        return VMCLocEnAndPoss<D, N, V>(wavef, poss, vps, useImpSamp, derivativeStep, masses, pot, coorBounds,
                                        numEnergies, gen);
    }};
    return VMCRBestParams_<D, N, V>(parBounds, wavef, enPossCalculator, numWalkers_gradDesc, function,
                                    boostrapSamples, gen);
}

//! @}

} // namespace vmcp

#endif
