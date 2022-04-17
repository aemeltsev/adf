// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#include "inc/complex.hpp"
namespace adf {

/**
 *  @brief - square root of complex number
 *  @param - unique pointer to complex value
 *  @return - returns the square root of complex number
 */
template<typename T>
inline complex<T> sqrt(complex<T>& p_val)
{
    T real = std::sqrt(p_val.mag()) * std::cos(p_val.arg()/2.);
    T imag = std::sqrt(p_val.mag()) * std::sin(p_val.arg()/2.);
    return complex<T>(real, imag);
}

/**
 *  @brief - factors quadratic equation with cmplx coefficients
 *  @param - a*x^2
 *  @param - b*x
 *  @param - c
 *  @return - std::pairs with complex roots
 */
template<typename T>
inline std::pair<complex<T>, complex<T>> quadr(
        complex<T>& p_a,
        complex<T>& p_b,
        complex<T>& p_c)
{
    complex<T> a2_var,
               ac4_tmp, ac4_var,
               sq1_tmp, sq2_tmp, sq_var,
               fr_root, sc_root,
               pb_tmp = p_b;
    complex<T> fr_const(2., 0.);
    complex<T> sc_const(4., 0.);

    a2_var = fr_const * p_a; /**< 2*a */

    ac4_tmp = sc_const * p_a;
    ac4_var = ac4_tmp * p_c; /**< 4*a*c */

    sq1_tmp = p_b * p_b;
    sq2_tmp = sq1_tmp - ac4_var;
    sq_var = adf::sqrt(sq2_tmp); /**< sqrt(b*b - 4*a*c) */

    pb_tmp = -pb_tmp;
    fr_root = pb_tmp + sq_var; /**< for first root */
    sc_root = pb_tmp - sq_var; /**< for second root */
    return std::make_pair(fr_root/a2_var, sc_root/a2_var);
}

template<typename T>
T abs(const complex<T> other)
{
    T out, re, im;
    re = std::abs(other.getReal());
    im = std::abs(other.getImag());
    out = std::min(re, im);
    re = std::max(re, im);
    out /= re;
    return re * std::sqrt(1. + out * out);
}

}
