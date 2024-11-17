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
        RandomGenerator gen{(std::random_device())()};
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
        vmcp::IntType const numberEnergies = 100;
        vmcp::CoordBounds<1> const bounds = {vmcp::Bound{vmcp::Coordinate{-100}, vmcp::Coordinate{100}}};
        vmcp::RandomGenerator rndGen{1};
        struct PotHO {
            vmcp::Mass m;
            vmcp::FPType omega;
            vmcp::FPType operator()(vmcp::Positions<1, 1> x) const {
                return x[0][0].val * x[0][0].val * (m.val * omega * omega / 2);
            }
        };
        Mass mass{1.f};
        FPType omega{2.6f};
        vmcp::VarParam bestAlpha{mass.val * omega / vmcp::hbar};
        vmcp::ParamBounds<1> alphaBound{
            vmcp::Bound{vmcp::VarParam{bestAlpha.val * 0.1f}, vmcp::VarParam{bestAlpha.val * 10}}};
        auto const wavefHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
            return std::exp(-alpha[0].val * x[0][0].val * x[0][0].val);
        }};
        auto const secondDerHO{[](vmcp::Positions<1, 1> x, vmcp::VarParams<1> alpha) {
            return (std::pow(x[0][0].val * alpha[0].val, 2) - alpha[0].val) *
                   std::exp(-alpha[0].val * x[0][0].val * x[0][0].val);
        }};
        PotHO potHO{mass, omega};
        vmcp::VMCResult const vmcr = vmcp::VMCEnergy<1, 1, 1>(wavefHO, alphaBound, secondDerHO, mass, potHO,
                                                              bounds, numberEnergies, rndGen);

        std::cout << vmcr.energy.val << '\t' << std::sqrt(vmcr.variance.val) << '\n';
    }
}
