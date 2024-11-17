//
//
// Contains the definition of the templates declared in vmcalgs.hpp
// This file is supposed to be #included at the end of vmcalgs.hpp and nowhere else
// It is just a way to improve the readability of vmcalgs.hpp
//

// FP/LF TODO: Couldn't just define everything in the header?

// Proposal of a rule: if it the function body is longer than one line, define it here

#ifndef VMCPROJECT_TYPES_INL
#define VMCPROJECT_TYPES_INL

#include "types.hpp"

namespace vmcp {

inline Coordinate &Coordinate::operator+=(Coordinate other) {
    val += other.val;
    return *this;
}
inline Coordinate &Coordinate::operator-=(Coordinate other) {
    val -= other.val;
    return *this;
}
inline Coordinate &Coordinate::operator*=(FPType other) {
    val *= other;
    return *this;
}
inline Coordinate &Coordinate::operator/=(FPType other) {
    val /= other;
    return *this;
}
inline Coordinate operator+(Coordinate lhs, Coordinate rhs) { return lhs += rhs; }
inline Coordinate operator-(Coordinate lhs, Coordinate rhs) { return lhs -= rhs; }
inline Coordinate operator*(Coordinate lhs, FPType rhs) { return lhs *= rhs; }
inline Coordinate operator/(Coordinate lhs, FPType rhs) { return lhs /= rhs; }

inline VarParam &VarParam::operator+=(VarParam other) {
    val += other.val;
    return *this;
}
inline VarParam &VarParam::operator-=(VarParam other) {
    val -= other.val;
    return *this;
}
inline VarParam &VarParam::operator*=(FPType other) {
    val *= other;
    return *this;
}
inline VarParam &VarParam::operator/=(FPType other) {
    val /= other;
    return *this;
}
inline VarParam operator+(VarParam lhs, VarParam rhs) { return lhs += rhs; }
inline VarParam operator-(VarParam lhs, VarParam rhs) { return lhs -= rhs; }
inline VarParam operator*(VarParam lhs, FPType rhs) { return lhs *= rhs; }
inline VarParam operator*(FPType lhs, VarParam rhs) { return rhs * lhs; }
inline VarParam operator/(VarParam lhs, FPType rhs) { return lhs /= rhs; }

inline bool operator<(Energy lhs, Energy rhs) { return lhs.val < rhs.val; }

} // namespace vmcp

#endif
