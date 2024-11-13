#include "src/vmcp.hpp"

#include <iomanip>
#include <iostream>
#include <numbers>

using namespace vmcp;

// The various features of main can be toggled here
bool const ho = true;
bool const statistics = true;
constexpr std::array features = {ho, statistics};

// LF TODO:
int main() {
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
        Bounds<1> const bounds = {Bound{-100, 100}};
        RandomGenerator gen{(std::random_device())()};
        VarParams<1> const initialAlpha{0.5f};
        Mass const mass{0.5f};

        for (FPType alphaVal = 0.1f; alphaVal <= 2; alphaVal += FPType{0.05f}) {
            VMCResult const vmcr = AvgAndVar_(Energies_(VMCEnAndPoss<1, 1, 1>(
                wavefHO, VarParams<1>{alphaVal}, secondDerHO, mass, potHO, bounds, numberEnergies, gen)));
            std::cout << "alpha: " << std::setprecision(3) << alphaVal << "\tenergy: " << std::setprecision(5)
                      << vmcr.energy.val << " +/- " << std::sqrt(vmcr.variance.val) << '\n';
        }
        VMCResult const vmcrBest =
            VMCEnergy<1, 1, 1>(wavefHO, initialAlpha, secondDerHO, mass, potHO, bounds, numberEnergies, gen);
        std::cout << "Energy with the best alpha:\n"
                  << std::setprecision(3) << "Energy: " << std::setprecision(5) << vmcrBest.energy.val
                  << " +/- " << std::sqrt(vmcrBest.variance.val) << '\n';

        if constexpr (features[1]) {
            FPType alphaVal = 0.9f;
            std::vector<Energy> energySamp = Energies_(VMCEnAndPoss<1, 1, 1>(
                wavefHO, VarParams<1>{alphaVal}, secondDerHO, mass, potHO, bounds, numberEnergies, gen));

            BlockingResult blockingResult = BlockingAnalysis(energySamp);
            for (size_t blockGroup = 0; blockGroup < blockingResult.sizes.size(); ++blockGroup) {
                std::cout << '\n'
                          << "block size: " << blockingResult.sizes[blockGroup]
                          << " , mean: " << blockingResult.means[blockGroup] << " , std. dev.: " << std::fixed
                          << std::setprecision(5) << blockingResult.stdDevs[blockGroup] << '\n';
            }

            UIntType const numSamples = 1000;
            BootstrapResult bootstrapResult = BootstrapAnalysis(energySamp, numSamples, gen);
            std::cout << '\n'
                      << "Bootstrap mean: " << bootstrapResult.meanOfMeans << std::fixed
                      << std::setprecision(5) << "\nBootstrap std. dev.: " << bootstrapResult.stdDevOfMeans
                      << "\nConfidence interval: " << bootstrapResult.confInterval.min << " - "
                      << bootstrapResult.confInterval.max << '\n';
        }
    }
}
