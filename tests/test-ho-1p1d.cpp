//!
//! @file test-ho-1p1d.cpp
//! @brief Tests for the harmonic oscillator system, with one particle in one dimension
//! @authors Lorenzo Fabbri, Francesco Orso Pancaldi
//!

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "test.hpp"
#include "vmcp.hpp"

#include <chrono>
#include <fstream>
#include <numbers>
#include <string>
#include <tuple>

TEST_CASE("Testing the harmonic oscillator") {
    std::ofstream file_stream;
    file_stream.open(logFilePath, std::ios_base::app);

    SUBCASE("One particle in one dimension") {
        vmcp::CoordBounds<1> const coordBound = {vmcp::Bound{vmcp::Coordinate{-100}, vmcp::Coordinate{100}}};
        vmcp::RandomGenerator rndGen{seed};
        vmcp::Masses<1> const mInit{1.f};
        vmcp::FPType const omegaInit = 1;
        vmcp::Mass const mStep{0.1f};
        vmcp::FPType const omegaStep = 0.1f;
        vmcp::IntType const mIterations = iterations;
        vmcp::IntType const omegaIterations = iterations;
        struct PotHO {
            vmcp::Mass m;
            vmcp::FPType omega;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
                return x[0][0].val * x[0][0].val * m.val * omega * omega / 2;
            }
        };
        vmcp::FPType const derivativeStep = coordBound[0].Length().val / derivativeStepDenom;

        SUBCASE("No variational parameters") {
            struct WavefHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return std::exp(-x[0][0].val * x[0][0].val * m.val * omega / (2 * vmcp::hbar));
                }
            };
            struct FirstDerHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return -x[0][0].val * m.val * omega / vmcp::hbar *
                           WavefHO{m, omega}(x, vmcp::VarParams<0>{});
                }
            };
            struct LaplHO {
                vmcp::Mass m;
                vmcp::FPType omega;
                vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<0>) const {
                    return (std::pow(x[0][0].val * m.val * omega / vmcp::hbar, 2) -
                            m.val * omega / vmcp::hbar) *
                           WavefHO{m, omega}(x, vmcp::VarParams<0>{});
                }
            };
            PotHO potHO{mInit[0], omegaInit};
            WavefHO wavefHO{mInit[0], omegaInit};
            vmcp::Gradients<1, 1, FirstDerHO> gradHO{mInit[0], omegaInit};
            vmcp::Laplacians<1, LaplHO> laplHO{mInit[0], omegaInit};

            auto start = std::chrono::high_resolution_clock::now();

            for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations; ++i, m_[0] += mStep) {
                potHO.m = m_[0];
                wavefHO.m = m_[0];
                gradHO[0][0].m = m_[0];
                laplHO[0].m = m_[0];
                for (auto [j, omega_] = std::tuple{vmcp::IntType{0}, omegaInit}; j != omegaIterations;
                     ++j, omega_ += omegaStep) {
                    potHO.omega = omega_;
                    wavefHO.omega = omega_;
                    gradHO[0][0].omega = omega_;
                    laplHO[0].omega = omega_;

                    vmcp::Energy const expectedEn{vmcp::hbar * omega_ / 2};
                    std::string const genericLogMes =
                        "mass: " + std::to_string(m_[0].val) + ", ang. vel.: " + std::to_string(omega_);
                    vmcp::Positions<1, 1> startPoss = FindPeak_<1, 1>(wavefHO, vmcp::VarParams<0>{}, potHO,
                                                                      coordBound, points_peakSearch, rndGen);

                    {
                        // Metropolis update, analytical derivative
                        std::string const logMes = metrLogMes + ", " + anDerLogMes + ", " + genericLogMes;
                        vmcp::VMCResult<0> const vmcr = vmcp::VMCEnergy<1, 1, 0>(
                            wavefHO, startPoss, vmcp::ParamBounds<0>{}, laplHO, m_, potHO, coordBound,
                            numEnergies, vmcp::StatFuncType::regular, bootstrapSamples, rndGen);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) < vmcEnergyTolerance, logMes);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) <
                                          max(vmcr.stdDev * allowedStdDevs, stdDevTolerance),
                                      logMes);
                    }
                    {
                        // Metropolis update, numerical derivative
                        std::string const logMes = metrLogMes + ", " + numDerLogMes + ", " + genericLogMes;
                        vmcp::VMCResult<0> const vmcr = vmcp::VMCEnergy<1, 1, 0>(
                            wavefHO, startPoss, vmcp::ParamBounds<0>{}, false, derivativeStep, m_, potHO,
                            coordBound, numEnergies, vmcp::StatFuncType::regular, bootstrapSamples, rndGen);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) < vmcEnergyTolerance, logMes);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) <
                                          max(vmcr.stdDev * allowedStdDevs, stdDevTolerance),
                                      logMes);
                    }
                    {
                        // Importance sampling update, analytical derivative
                        std::string const logMes = impSampLogMes + ", " + anDerLogMes + ", " + genericLogMes;
                        vmcp::VMCResult<0> const vmcr = vmcp::VMCEnergy<1, 1, 0>(
                            wavefHO, startPoss, vmcp::ParamBounds<0>{}, gradHO, laplHO, m_, potHO, coordBound,
                            numEnergies, vmcp::StatFuncType::regular, bootstrapSamples, rndGen);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) < vmcEnergyTolerance, logMes);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) <
                                          max(vmcr.stdDev * allowedStdDevs, stdDevTolerance),
                                      logMes);
                    }
                    /* {
                        // Importance sampling update, numerical derivative
                        std::string const logMes = impSampLogMes + ", " + numDerLogMes + ", " + genericLogMes;
                        vmcp::VMCResult<0> const vmcr = vmcp::VMCEnergy<1, 1, 0>(
                            wavefHO, startPoss, vmcp::ParamBounds<0>{}, true, derivativeStep, m_, potHO,
                    coordBound, numEnergies, vmcp::StatFuncType::regular, bootstrapSamples, rndGen);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) < vmcEnergyTolerance, logMes);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) <
                                          max(vmcr.stdDev * allowedStdDevs, stdDevTolerance),
                                      logMes);
                    } */
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::seconds>(stop - start);
            file_stream << "1p1d harmonic oscillator, no var. parameters (seconds): " << duration.count()
                        << '\n';
        }

        SUBCASE("One variational parameter") {
            PotHO potHO{mInit[0], omegaInit};
            auto const wavefHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
                return std::exp(-alpha[0].val * x[0][0].val * x[0][0].val / 2);
            }};
            std::array const laplHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
                return (std::pow(x[0][0].val * alpha[0].val, 2) - alpha[0].val) *
                       std::exp(-alpha[0].val * x[0][0].val * x[0][0].val / 2);
            }};

            auto start = std::chrono::high_resolution_clock::now();

            for (auto [i, m_] = std::tuple{vmcp::IntType{0}, mInit}; i != mIterations;
                 i += vpIterationsFactor, m_[0] += mStep * vpIterationsFactor) {
                potHO.m = m_[0];
                for (auto [j, omega_] = std::tuple{vmcp::IntType{0}, omegaInit}; j != omegaIterations;
                     j += vpIterationsFactor / 2, omega_ += omegaStep * vpIterationsFactor / 2) {
                    potHO.omega = omega_;

                    vmcp::VarParam const bestParam{m_[0].val * omega_ / vmcp::hbar};
                    vmcp::ParamBounds<1> const parBound{
                        NiceBound(bestParam, vmcp::minParamFactor, vmcp::maxParamFactor, vmcp::maxParDiff)};
                    vmcp::Energy const expectedEn{vmcp::hbar * omega_ / 2};
                    std::string const genericLogMes =
                        "mass: " + std::to_string(m_[0].val) + ", ang. vel.: " + std::to_string(omega_);
                    vmcp::Positions<1, 1> const startPoss = FindPeak_<1, 1>(
                        wavefHO, vmcp::VarParams<1>{bestParam}, potHO, coordBound, points_peakSearch, rndGen);

                    {
                        // Metropolis algorithm, analytical derivative
                        std::string const logMes = metrLogMes + ", " + anDerLogMes + ", " + genericLogMes;
                        auto startVP = std::chrono::high_resolution_clock::now();
                        vmcp::VMCResult<1> const vmcr =
                            vmcp::VMCEnergy<1, 1, 1>(wavefHO, startPoss, parBound, laplHO, std::array{m_},
                                                     potHO, coordBound, numEnergies / vpNumEnergiesFactor,
                                                     vmcp::StatFuncType::regular, bootstrapSamples, rndGen);
                        auto stopVP = std::chrono::high_resolution_clock::now();
                        auto durationVP = duration_cast<std::chrono::seconds>(stopVP - startVP);
                        file_stream << logMes << " (seconds): " << durationVP.count() << '\n';

                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) < vmcEnergyTolerance, logMes);
                        CHECK_MESSAGE(abs(vmcr.energy - expectedEn) <
                                          max(vmcr.stdDev * allowedStdDevs, stdDevTolerance),
                                      logMes);
                    }
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::seconds>(stop - start);
            file_stream << "1p1d harmonic oscillator, one var. parameter (seconds): " << duration.count()
                        << '\n';
        }
    }
}
