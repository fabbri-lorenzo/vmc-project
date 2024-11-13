#include "vmcalgs.inl"

namespace vmcp {

VMCResult AvgAndVar_(std::vector<Energy> const &v) {
    assert(v.size() > 1);
    auto const size = v.size();
    Energy const avg{std::accumulate(v.begin(), v.end(), Energy{0},
                                     [](Energy e1, Energy e2) { return Energy{e1.val + e2.val}; })
                         .val /
                     size};
    EnVariance const var{std::accumulate(v.begin(), v.end(), EnVariance{0},
                                         [avg](EnVariance ev, Energy e) {
                                             return EnVariance{ev.val + std::pow(e.val - avg.val, 2)};
                                         })
                             .val /
                         (size - 1)};
    return VMCResult{avg, var};
}

} // namespace vmcp
