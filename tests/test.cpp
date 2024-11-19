#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "vmcp.hpp"

#include <chrono>
#include <fstream>
#include <functional>
#include <numbers>
#include <string>
#include <tuple>

// LF TODO: Rename the gradients
// I will explain myself better

// FP TODO: Change mStep to a vmcp::Mass everywhere and define operator+= for vmcp::Mass

// Chosen at random
constexpr vmcp::UIntType seed = 648263u;
constexpr vmcp::IntType iterations = 32;
// FP TODO: Rename this
constexpr vmcp::IntType allowedStdDevs = 5;
// FP TODO: Explain, and maybe rename
constexpr vmcp::FPType stdDevTolerance = 1e-9f;
// FP TODO: Rename
constexpr vmcp::IntType varParamsFactor = 16;
// FP TODO: One pair of brackets can probably be removed here
// Learn the priority of the operations and adjust the rest of the code too
static_assert((iterations % varParamsFactor) == 0);
// LF TODO: This is unused for now!
constexpr bool testImpSamp = true;
// Rules out situations where both the VMC energy and the variance are extremely large
constexpr vmcp::FPType vmcEnergyTolerance = 0.5f;
// FP TODO: Explain
constexpr vmcp::FPType minParamFactor = 0.33f;
constexpr vmcp::FPType maxParamFactor = 3;
const std::string logFileName = "../artifacts/test-log.txt";

TEST_CASE("Testing VMCLocEnAndPoss_") {
    std::ofstream file_stream;
    file_stream.open(logFileName, std::ios_base::app);

    SUBCASE("1D harmonic oscillator") {
        vmcp::IntType const numberEnergies = 100;
        vmcp::CoordBounds<1> const coordBound = {vmcp::Bound{vmcp::Coordinate{-100}, vmcp::Coordinate{100}}};
        vmcp::RandomGenerator rndGen{seed};
        vmcp::Mass const mInit{1.f};
        vmcp::FPType const omegaInit = 1.f;
        vmcp::FPType const mStep = 0.2f;
        vmcp::FPType const omegaStep = 0.2f;
        vmcp::IntType const mIterations = iterations;
        vmcp::IntType const omegaIterations = iterations;
        struct PotHO {
            vmcp::Mass m;
            vmcp::FPType omega;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
                return x[0][0].val * x[0][0].val * (m.val * omega * omega / 2);
            }
        };

        SUBCASE("No variational parameters, with Metropolis or importance sampling") {
            struct WavefHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return std::exp(-x[0][0].val * x[0][0].val * (m.val * omega / (2 * vmcp::hbar)));
                }
            };
            struct FirstDerHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return (-x[0][0].val * (m.val * omega / (vmcp::hbar))) *
                           WavefHO{m, omega}(x, vmcp::VarParams<0>{});
                }
            };
            struct LaplHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return (std::pow(x[0][0].val * m.val * omega / vmcp::hbar, 2) -
                            (m.val * omega / vmcp::hbar)) *
                           WavefHO{m, omega}(x, vmcp::VarParams<0>{});
                }
            };

            WavefHO wavefHO{mInit.val, omegaInit};
            vmcp::Gradients<1, 1, FirstDerHO> gradHO{FirstDerHO{mInit.val, omegaInit}};
            vmcp::Laplacians<1, LaplHO> laplHO{LaplHO{mInit, omegaInit}};
            PotHO potHO{mInit.val, omegaInit};

            auto start = std::chrono::high_resolution_clock::now();

            for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations; ++i, m_.val += mStep) {
                potHO.m = m_;
                wavefHO.m = m_;
                gradHO[0][0].m = m_;
                laplHO[0].m = m_;
                for (auto [j, omega_] = std::tuple{vmcp::IntType{0}, omegaInit}; j != omegaIterations;
                     ++j, omega_ += omegaStep) {
                    wavefHO.omega = omega_;
                    potHO.omega = omega_;
                    laplHO[0].omega = omega_;

                    vmcp::Energy const expectedEn{vmcp::hbar * omega_ / 2};
                    vmcp::VMCResult const vmcrMetr =
                        vmcp::VMCEnergy<1, 1, 0>(wavefHO, vmcp::ParamBounds<0>{}, laplHO, m_, potHO,
                                                 coordBound, numberEnergies, rndGen);
                    vmcp::VMCResult const vmcrImpSamp =
                        vmcp::VMCEnergy<1, 1, 0>(wavefHO, vmcp::ParamBounds<0>{}, gradHO, laplHO, m_, potHO,
                                                 coordBound, numberEnergies, rndGen);

                    std::string logMessage{"mass: " + std::to_string(m_.val) +
                                           ", ang. vel.: " + std::to_string(omega_)};
                    CHECK_MESSAGE(std::abs(vmcrMetr.energy.val - expectedEn.val) < vmcEnergyTolerance,
                                  logMessage);
                    CHECK_MESSAGE(
                        std::abs(vmcrMetr.energy.val - expectedEn.val) <
                            std::max((allowedStdDevs * std::sqrt(vmcrMetr.variance.val)), stdDevTolerance),
                        logMessage);
                    CHECK_MESSAGE(std::abs(vmcrImpSamp.energy.val - expectedEn.val) < vmcEnergyTolerance,
                                  logMessage);
                    CHECK_MESSAGE(
                        std::abs(vmcrImpSamp.energy.val - expectedEn.val) <
                            std::max((allowedStdDevs * std::sqrt(vmcrImpSamp.variance.val)), stdDevTolerance),
                        logMessage);
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::seconds>(stop - start);
            file_stream << "Harmonic oscillator, no var. parameters (seconds): " << duration.count() << '\n';
        }

        SUBCASE("One variational parameter") {
            PotHO potHO{mInit, omegaInit};
            auto const wavefHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
                return std::exp(-alpha[0].val * x[0][0].val * x[0][0].val);
            }};
            std::array laplHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
                return (std::pow(x[0][0].val * alpha[0].val, 2) - alpha[0].val) *
                       std::exp(-alpha[0].val * x[0][0].val * x[0][0].val);
            }};

            auto start = std::chrono::high_resolution_clock::now();

            for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations;
                 i += varParamsFactor, m_.val += mStep * varParamsFactor) {
                potHO.m = m_;
                for (auto [j, omega_] = std::tuple{vmcp::IntType{0}, omegaInit}; j != omegaIterations;
                     j += varParamsFactor, omega_ += omegaStep * varParamsFactor) {
                    potHO.omega = omega_;

                    vmcp::VarParam bestParam{m_.val * omega_ / vmcp::hbar};
                    vmcp::ParamBounds<1> const parBound{
                        vmcp::Bound{bestParam * minParamFactor, bestParam * maxParamFactor}};
                    vmcp::Energy const expectedEn{vmcp::hbar * omega_ / 2};

                    auto startOnePar = std::chrono::high_resolution_clock::now();
                    vmcp::VMCResult const vmcr = vmcp::VMCEnergy<1, 1, 1>(
                        wavefHO, parBound, laplHO, m_, potHO, coordBound, numberEnergies, rndGen);
                    auto stopOnePar = std::chrono::high_resolution_clock::now();
                    auto durationOnePar = duration_cast<std::chrono::seconds>(stopOnePar - startOnePar);
                    file_stream << "Harmonic oscillator, one var. parameter, with mass " << m_.val
                                << " and ang. vel. " << omega_ << " (seconds): " << durationOnePar.count()
                                << '\n';

                    std::string logMessage{"mass: " + std::to_string(m_.val) +
                                           ", ang. vel.: " + std::to_string(omega_)};
                    CHECK_MESSAGE(std::abs(vmcr.energy.val - expectedEn.val) < vmcEnergyTolerance,
                                  logMessage);
                    CHECK_MESSAGE(
                        std::abs(vmcr.energy.val - expectedEn.val) <
                            std::max((allowedStdDevs * std::sqrt(vmcr.variance.val)), stdDevTolerance),
                        logMessage);
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::seconds>(stop - start);
            file_stream << "Harmonic oscillator, one var. parameter (seconds): " << duration.count() << '\n';
        }
    }

    SUBCASE("1D potential box") {
        // l = length of the box
        vmcp::IntType const numberEnergies = 100;
        vmcp::RandomGenerator rndGen{seed};
        vmcp::Mass const mInit{1.f};
        vmcp::FPType const lInit = 1;
        vmcp::FPType const mStep = 0.2f;
        vmcp::FPType const lStep = 0.2f;
        vmcp::IntType const mIterations = iterations;
        vmcp::IntType const lIterations = iterations;
        // Use the following potential: V(x)
        // = 0 if |x| < 99/100 * l/2
        // = V_0 * (200/l)^2 * (|x| - 99/100 * l/2) if |x| >= 99/100 * l/2
        struct PotBox {
            vmcp::FPType l;
            vmcp::FPType V_0;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
                if (std::abs(x[0][0].val < 99 * l / 200)) {
                    return 0;
                } else {
                    return std::pow((20 / l) * (std::abs(x[0][0].val) - 9 * l / 20), 2);
                }
            }
        };
        // FP TODO: Fix this (related to the 1 var par case)
        // struct PotBox {
        //    vmcp::FPType l;
        //    vmcp::FPType V_0;
        //    vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
        //        return 2 * V_0 / (1 + std::exp(-(20 * std::log(9) / l) * (std::abs(x[0][0].val) - l / 2)));
        //    }
        // };

        SUBCASE("No variational parameters, with Metropolis or importance sampling") {
            struct WavefBox {
                vmcp::FPType l;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    if (std::abs(x[0][0].val) <= l / 2) {
                        return std::cos(std::numbers::pi_v<vmcp::FPType> * x[0][0].val / l);
                    } else {
                        return 0;
                    }
                }
            };
            struct FirstDerBox {
                vmcp::FPType l;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    if (std::abs(x[0][0].val) <= l / 2) {
                        return -(std::numbers::pi_v<vmcp::FPType> / l) *
                               std::sin(std::numbers::pi_v<vmcp::FPType> * x[0][0].val / l);
                    } else {
                        return 0;
                    }
                }
            };
            struct LaplBox {
                vmcp::FPType l;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    if (std::abs(x[0][0].val) <= l / 2) {
                        return -std::pow(std::numbers::pi_v<vmcp::FPType> / l, 2) *
                               std::cos(std::numbers::pi_v<vmcp::FPType> * x[0][0].val / l);
                    } else {
                        return 0;
                    }
                }
            };

            WavefBox wavefBox{lInit};
            vmcp::Gradients<1, 1, FirstDerBox> gradBox{FirstDerBox{lInit}};
            vmcp::Laplacians<1, LaplBox> laplBox{lInit};
            PotBox potBox{0, lInit};

            auto start = std::chrono::high_resolution_clock::now();

            for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations; ++i, m_.val += mStep) {
                for (auto [j, l_] = std::tuple{vmcp::IntType{0}, lInit}; j != lIterations; ++j, l_ += lStep) {
                    potBox.l = l_;
                    wavefBox.l = l_;
                    gradBox[0][0].l = l_;
                    laplBox[0].l = l_;

                    vmcp::CoordBounds<1> const coorBound{
                        vmcp::Bound{vmcp::Coordinate{-l_ / 2}, vmcp::Coordinate{l_ / 2}}};
                    vmcp::Energy const expectedEn{
                        1 / (2 * m_.val) * std::pow(vmcp::hbar * std::numbers::pi_v<vmcp::FPType> / l_, 2)};
                    potBox.V_0 = 10 * expectedEn.val;
                    vmcp::VMCResult const vmcrMetr =
                        vmcp::VMCEnergy<1, 1, 0>(wavefBox, vmcp::ParamBounds<0>{}, laplBox, m_, potBox,
                                                 coorBound, numberEnergies, rndGen);
                    vmcp::VMCResult const vmcrImpSamp =
                        vmcp::VMCEnergy<1, 1, 0>(wavefBox, vmcp::ParamBounds<0>{}, gradBox, laplBox, m_,
                                                 potBox, coorBound, numberEnergies, rndGen);

                    std::string logMessage{"mass: " + std::to_string(m_.val) +
                                           ", length: " + std::to_string(l_)};
                    CHECK_MESSAGE(std::abs(vmcrMetr.energy.val - expectedEn.val) < vmcEnergyTolerance,
                                  logMessage);
                    CHECK_MESSAGE(
                        std::abs(vmcrMetr.energy.val - expectedEn.val) <
                            std::max((allowedStdDevs * std::sqrt(vmcrMetr.variance.val)), stdDevTolerance),
                        logMessage);
                    CHECK_MESSAGE(std::abs(vmcrImpSamp.energy.val - expectedEn.val) < vmcEnergyTolerance,
                                  logMessage);
                    CHECK_MESSAGE(
                        std::abs(vmcrImpSamp.energy.val - expectedEn.val) <
                            std::max((allowedStdDevs * std::sqrt(vmcrImpSamp.variance.val)), stdDevTolerance),
                        logMessage);
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::seconds>(stop - start);
            file_stream << "Particle in a box, no var. parameters (seconds): " << duration.count() << '\n';
        }

        // FP TODO: Fix this
        // SUBCASE("One variational parameter") {
        //     struct WavefBox {
        //         vmcp::FPType l;
        //         vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) const {
        //             if (std::abs(x[0][0].val) <= l / 2) {
        //                 return std::abs(std::cos(alpha[0].val * x[0][0].val)) +
        //                        0.001 * std::abs(std::cos(alpha[0].val * l / 2)) /
        //                            (vmcp::FPType{10.f} - std::pow(2 * x[0][0].val / l, 2));
        //             } else {
        //                 return 0;
        //             }
        //         }
        //     };
        //     struct SecondDerBox {
        //         vmcp::FPType l;
        //         vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) const {
        //             return -alpha[0].val * alpha[0].val * WavefBox{l}(x, alpha);
        //         }
        //     };
        //
        //    WavefBox wavefBox{lInit};
        //    SecondDerBox secondDerBox{lInit};
        //    PotBox potBox{0, lInit};
        //
        //    auto start = std::chrono::high_resolution_clock::now();
        //
        //    for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations;
        //         i += varParamsFactor, m_.val += mStep * varParamsFactor) {
        //        for (auto [j, l_] = std::tuple{vmcp::IntType{0}, lInit}; j != lIterations;
        //             j += varParamsFactor, l_ += lStep * varParamsFactor) {
        //            wavefBox.l = l_;
        //            secondDerBox.l = l_;
        //            potBox.l = l_;
        //
        //            vmcp::CoordBounds<1> const coordBound{
        //                vmcp::Bound{vmcp::Coordinate{-l_ / 2}, vmcp::Coordinate{l_ / 2}}};
        //            vmcp::VarParam bestParam{std::numbers::pi_v<vmcp::FPType> / l_};
        //            vmcp::ParamBounds<1> const parBound{
        //                vmcp::Bound{bestParam * minParamFactor, bestParam * maxParamFactor}};
        //
        //            vmcp::Energy expectedEn{1 / (2 * m_.val) *
        //                                    std::pow(vmcp::hbar * std::numbers::pi_v<vmcp::FPType> / l_,
        //                                    2)};
        //            potBox.V_0 = 5 * expectedEn.val;
        //
        //            auto startOnePar = std::chrono::high_resolution_clock::now();
        //            vmcp::VMCResult const vmcr = vmcp::VMCEnergy<1, 1, 1>(
        //                wavefBox, parBound, secondDerBox, m_, potBox, coordBound, numberEnergies, rndGen);
        //            auto stopOnePar = std::chrono::high_resolution_clock::now();
        //            auto durationOnePar = duration_cast<std::chrono::seconds>(stopOnePar - startOnePar);
        //            file_stream << "Harmonic oscillator, one var. parameter, with mass " << m_.val
        //                        << " and length " << l_ << " (seconds): " << durationOnePar.count() << '\n';
        //
        //            CHECK(std::abs(vmcr.energy.val - expectedEn.val) < vmcEnergyTolerance);
        //            CHECK(std::abs(vmcr.energy.val - expectedEn.val) <
        //                  std::max(allowedStdDevs * std::sqrt(vmcr.variance.val), stdDevTolerance));
        //        }
        //    }
        //
        //    auto stop = std::chrono::high_resolution_clock::now();
        //    auto duration = duration_cast<std::chrono::seconds>(stop - start);
        //    file_stream << "Particle in a box, one var. parameter (seconds): " << duration.count() << '\n';
        //}
    }

    // LF TODO: to be finished & adding references
    // SUBCASE("1D triangular well potential") {
    //     struct WavefTrg {
    //         vmcp::FPType alpha;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
    //             return std::pow(std::numbers::e_v<vmcp::FPType>,
    //                             -std::pow(2 * m * alpha / (vmcp::hbar * vmcp::hbar), 1 / 3) *
    //                                 std::abs(x[0][0].val));
    //         }
    //     };
    //     struct PotTrg {
    //         vmcp::FPType alpha;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> x) const { return alpha *
    //         std::abs(x[0][0].val);
    //         }
    //     };
    //     struct GradTrg {
    //         vmcp::FPType alpha;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
    //             return -std::pow(2 * m * alpha / (vmcp::hbar * vmcp::hbar), 1 / 3) *
    //                    std::signbit(x[0][0].val) * WavefTrg{alpha, m}(x, vmcp::VarParams<0>{});
    //         }
    //     };
    //     struct SecondDerTrg {
    //         vmcp::FPType alpha;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
    //             return std::pow(2 * m * alpha / (vmcp::hbar * vmcp::hbar), 2 / 3) *
    //                    WavefTrg{alpha, m}(x, vmcp::VarParams<0>{});
    //         }
    //     };
    //
    //    vmcp::IntType const numberEnergies = 100;
    //    vmcp::Bounds<1> const bounds = {vmcp::Bound{-100, 100}};
    //    vmcp::RandomGenerator rndGen{seed};
    //    vmcp::FPType const alphaInit = 0.1f;
    //    vmcp::FPType const mInit = 0.5f;
    //    vmcp::FPType const alphaStep = 0.2f;
    //    vmcp::FPType const mStep = 0.2f;
    //    vmcp::IntType const alphaIterations = 10;
    //    vmcp::IntType const mIterations = 10;
    //    vmcp::Mass mass = vmcp::Mass{1.f};
    //    WavefTrg wavefTrg{alphaInit, mInit};
    //    PotTrg potTrg{alphaInit};
    //    GradTrg gradTrg{alphaInit, mInit};
    //    SecondDerTrg secondDerTrg{alphaInit, mInit};
    //    std::array<GradTrg, 1> gradTrgArr = {gradTrg};
    //    for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations; ++i, m_ += mStep) {
    //        wavefTrg.m = m_;
    //        gradTrg.m = m_;
    //        secondDerTrg.m = m_;
    //        for (auto [j, alpha_] = std::tuple{vmcp::IntType{0}, alphaInit}; j != alphaIterations;
    //             ++j, alpha_ += alphaStep) {
    //            wavefTrg.alpha = alpha_;
    //            potTrg.alpha = alpha_;
    //            gradTrg.alpha = alpha_;
    //            secondDerTrg.alpha = alpha_;
    //            mass.val = m_;
    //
    //            vmcp::Energy expectedEn{std::pow(2 * m_, -1 / 3) *
    //                                    std::pow(alpha_ * vmcp::hbar * 3 / 4, 2 / 3)};
    //            vmcp::VMCResult vmcr =
    //                vmcp::VMCEnergy<1, 1>(wavefTrg, vmcp::VarParams<0>{}, testImpSamp, gradTrgArr,
    //                                      secondDerTrg, mass, potTrg, bounds, numberEnergies, rndGen);
    //
    //            if (std::abs(vmcr.variance.val) < varianceTolerance) {
    //                CHECK(std::abs(vmcr.energy.val - expectedEn.val) < varianceTolerance);
    //            } else {
    //                CHECK(std::abs(vmcr.energy.val - expectedEn.val) < vmcr.variance.val);
    //            }
    //        }
    //    }
    //}

    // LF TODO: to be finished & adding references
    // SUBCASE("1D radial Schroedinger eq") {
    //     struct WavefRad {
    //         vmcp::FPType k;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> r, vmcp::VarParams<0>) const {
    //             return r[0][0].val * std::pow(std::numbers::e_v<vmcp::FPType>,
    //                                           -m * k / (vmcp::hbar * vmcp::hbar) * r[0][0].val);
    //         }
    //     };
    //     struct PotRad {
    //         vmcp::FPType k;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> r) const { return (-k / r[0][0].val); }
    //     };
    //     struct GradRad {
    //         vmcp::FPType k;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> r, vmcp::VarParams<0>) const {
    //             return (1 - m * k / (vmcp::hbar * vmcp::hbar) * r[0][0].val) *
    //                    WavefRad{k, m}(r, vmcp::VarParams<0>{});
    //         }
    //     };
    //     struct SecondDerRad {
    //         vmcp::FPType k;
    //         vmcp::FPType m;
    //         vmcp::FPType operator()(vmcp::Positions<1, 1> r, vmcp::VarParams<0>) const {
    //             return (std::pow(1 - r[0][0].val * m * k / (vmcp::hbar * vmcp::hbar), 2) -
    //                     m * k / (vmcp::hbar * vmcp::hbar)) *
    //                    WavefRad{k, m}(r, vmcp::VarParams<0>{});
    //         }
    //     };
    //
    //    vmcp::IntType const numberEnergies = 100;
    //    vmcp::Bounds<1> const bounds = {vmcp::Bound{0, 100}};
    //    vmcp::RandomGenerator rndGen{seed};
    //    vmcp::FPType const kInit = 1.2f;
    //    vmcp::FPType const mInit = 0.3f;
    //    vmcp::FPType const kStep = 0.1f;
    //    vmcp::FPType const mStep = 0.2f;
    //    vmcp::IntType const kIterations = 10;
    //    vmcp::IntType const mIterations = 10;
    //    vmcp::Mass mass = vmcp::Mass{1.f};
    //    WavefRad wavefRad{kInit, mInit};
    //    PotRad potRad{kInit};
    //    GradRad gradRad{kInit, mInit};
    //    SecondDerRad secondDerRad{kInit, mInit};
    //    std::array<GradRad, 1> gradRadArr = {gradRad};
    //
    //    for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations; ++i, m_ += mStep) {
    //        wavefRad.m = m_;
    //        gradRad.m = m_;
    //        secondDerRad.m = m_;
    //        for (auto [j, k_] = std::tuple{vmcp::IntType{0}, kInit}; j != kIterations; ++j, k_ += kStep)
    //        {
    //            wavefRad.k = k_;
    //            potRad.k = k_;
    //            gradRad.k = k_;
    //            secondDerRad.k = k_;
    //            mass.val = m_;
    //
    //            vmcp::Energy expectedEn{-m_ * std::pow((1.6021766 * std::pow(10, -19)), 4) /
    //                                    (2 * vmcp::hbar)};
    //            vmcp::VMCResult vmcr =
    //                vmcp::VMCEnergy<1, 1>(wavefRad, vmcp::VarParams<0>{}, testImpSamp, gradRadArr,
    //                                      secondDerRad, mass, potRad, bounds, numberEnergies, rndGen);
    //
    //            if (std::abs(vmcr.variance.val) < varianceTolerance) {
    //                CHECK(std::abs(vmcr.energy.val - expectedEn.val) < varianceTolerance);
    //            } else {
    //                CHECK(std::abs(vmcr.energy.val - expectedEn.val) < vmcr.variance.val);
    //            }
    //        }
    //    }
    //}
}
