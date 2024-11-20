#include "src/vmcp.hpp"

#include <iomanip>
#include <iostream>
#include <numbers>

using namespace vmcp;

// The various features of main can be toggled here
bool const ho = true;
bool const statistics = true;
bool const bugfixing = false;
constexpr std::array features = {ho, statistics, bugfixing};


int main() {
    RandomGenerator gen{(std::random_device())()};

    // Feature 1:
    // Various tests with the harmonic oscillator
    if constexpr (features[0]) {
        // It is not normalized, but it doesn't matter
        auto const wavefHO{[](Positions<1, 1> x, VarParams<1> alpha) {
            return std::pow(std::numbers::e_v<FPType>, -alpha[0].val * x[0][0].val * x[0][0].val / 2);
        }};
        auto const potHO{[](Positions<1, 1> x) { return x[0][0].val * x[0][0].val; }};
        auto const secondDerHO{[&wavefHO](Positions<1, 1> x, VarParams<1> alpha) {
            return (std::pow(alpha[0].val * x[0][0].val, 2) - alpha[0].val) * wavefHO(x, alpha);
        }};

        int const numberEnergies = 100;
        CoordBounds<1> const coorBounds = {Bound{Coordinate{-100}, Coordinate{100}}};
        Mass const mass{0.5f};

        for (VarParam alpha{0.1f}; alpha.val <= 2; alpha.val += FPType{0.05f}) {
            VMCResult const vmcr = AvgAndVar_(LocalEnergies_(VMCLocEnAndPoss<1, 1, 1>(
                wavefHO, VarParams<1>{alpha}, secondDerHO, mass, potHO, coorBounds, numberEnergies, gen)));
            std::cout << "alpha: " << std::setprecision(3) << alpha.val
                      << "\tenergy: " << std::setprecision(5) << vmcr.energy.val << " +/- "
                      << std::sqrt(vmcr.variance.val) << '\n';
        }

        ParamBounds<1> alphaBounds{Bound{VarParam{0.5f}, VarParam{1.5f}}};

        VMCResult const vmcrBest = VMCEnergy<1, 1, 1>(wavefHO, alphaBounds, secondDerHO, mass, potHO,
                                                      coorBounds, numberEnergies, gen);
        std::cout << "Energy with the best alpha:\n"
                  << std::setprecision(3) << "Energy: " << std::setprecision(5) << vmcrBest.energy.val
                  << " +/- " << std::sqrt(vmcrBest.variance.val) << '\n';
    }
    // Feature 2
    // Statistical analysis
    if constexpr (features[1]) {
        auto const wavefHO{[](Positions<1, 1> x, VarParams<1> alpha) {
            return std::pow(std::numbers::e_v<FPType>, -alpha[0].val * x[0][0].val * x[0][0].val / 2);
        }};
        auto const potHO{[](Positions<1, 1> x) { return x[0][0].val * x[0][0].val; }};
        auto const secondDerHO{[&wavefHO](Positions<1, 1> x, VarParams<1> alpha) {
            return (std::pow(alpha[0].val * x[0][0].val, 2) - alpha[0].val) * wavefHO(x, alpha);
        }};

        int const numberEnergies = 1000;
        CoordBounds<1> const coorBounds = {Bound{Coordinate{-100}, Coordinate{100}}};
        RandomGenerator gen{(std::random_device())()};
        Mass const mass{0.5f};
        FPType alphaVal = 0.9f;
        std::vector<Energy> energySamp = LocalEnergies_(VMCLocEnAndPoss<1, 1, 1>(
            wavefHO, VarParams<1>{alphaVal}, secondDerHO, mass, potHO, coorBounds, numberEnergies, gen));

        BlockingResult blockingResult = BlockingAnalysis(energySamp);
        for (size_t blockGroup = 0; blockGroup < blockingResult.sizes.size(); ++blockGroup) {
            std::cout << '\n'
                      << "block size: " << blockingResult.sizes[blockGroup]
                      << " , mean: " << blockingResult.means[blockGroup] << " , std. dev.: " << std::fixed
                      << std::setprecision(5) << blockingResult.stdDevs[blockGroup] << '\n';
        }

        UIntType const numSamples = 10000;
        BootstrapResult bootstrapResult = BootstrapAnalysis(energySamp, numSamples, gen);
        std::cout << '\n'
                  << "Bootstrap mean: " << bootstrapResult.mean << std::fixed << std::setprecision(5)
                  << "\nBootstrap std. dev.: " << bootstrapResult.stdDev
                  << "\nConfidence interval with confidence level of 95% : "
                  << bootstrapResult.confInterval.min << " - " << bootstrapResult.confInterval.max << '\n';
    }

    // Feature 3
    // Just bugfixing
    if constexpr (features[2]) {
        vmcp::Mass const mass{1.f};
        vmcp::FPType const length = 1;
        struct WavefBox {
            vmcp::FPType l;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) const {
                if (std::abs(x[0][0].val) <= l / 2) {
                    return std::abs(std::cos(alpha[0].val * x[0][0].val)) +
                           10 * std::abs(std::cos(alpha[0].val * l / 2)) /
                               (5 - std::pow(2 * x[0][0].val / l, 2));
                } else {
                    return 0;
                }
            }
        };
        struct PotBox {
            vmcp::FPType l;
            vmcp::FPType V_0;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
                return 2 * V_0 / (1 + std::exp(-(20 * std::log(9) / l) * (std::abs(x[0][0].val) - l / 2)));
            }
        };
        struct SecondDerBox {
            vmcp::FPType l;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) const {
                return -alpha[0].val * alpha[0].val * WavefBox{l}(x, alpha);
            }
        };

        vmcp::IntType const numberEnergies = 100;
        WavefBox wavefBox{length};
        SecondDerBox secondDerBox{length};
        vmcp::Energy const expectedEn{1 / (2 * mass.val) *
                                      std::pow(vmcp::hbar * std::numbers::pi_v<vmcp::FPType> / length, 2)};
        PotBox potBox{20 * expectedEn.val, length};
        vmcp::VarParam bestParam{std::numbers::pi_v<vmcp::FPType> / length};
        vmcp::CoordBounds<1> const coordBound{
            vmcp::Bound{vmcp::Coordinate{-length / 2}, vmcp::Coordinate{length / 2}}};

        for (vmcp::VarParam alpha = bestParam / 100; alpha.val <= 2 * bestParam.val;
             alpha.val += bestParam.val / 100) {
            VMCResult const vmcr = AvgAndVar_(LocalEnergies_(VMCLocEnAndPoss<1, 1, 1>(
                wavefBox, VarParams<1>{alpha}, secondDerBox, mass, potBox, coordBound, numberEnergies, gen)));
            std::cout << "alpha: " << std::setprecision(3) << alpha.val
                      << "\tenergy: " << std::setprecision(5) << vmcr.energy.val << " +/- "
                      << std::sqrt(vmcr.variance.val) << '\n';
        }
    }
}
