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
inline Coordinate operator+(Coordinate lhs, Coordinate rhs) {
    Coordinate result = lhs;
    return result += rhs;
}
inline Coordinate &Coordinate::operator-=(Coordinate other) {
    val -= other.val;
    return *this;
}
inline Coordinate operator-(Coordinate lhs, Coordinate rhs) {
    Coordinate result = lhs;
    return result -= rhs;
}


inline VarParam &VarParam::operator+=(VarParam other) {
    val += other.val;
    return *this;
}
inline VarParam operator+(VarParam lhs, VarParam rhs) {
    VarParam result = lhs;
    return result += rhs;
}
inline VarParam &VarParam::operator-=(VarParam other) {
    val -= other.val;
    return *this;
}
inline VarParam operator-(VarParam lhs, VarParam rhs) {
    VarParam result = lhs;
    return result -= rhs;
}




} // namespace vmcp

#endif
