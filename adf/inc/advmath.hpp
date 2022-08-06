// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef ADVMATH_H
#define ADVMATH_H
#include <cstdint>
#include <cmath>
#include <array>

#include "base.hpp"

namespace adf {
//! For more see Sidhu S.S.-Elliptic integrals and functions.- Computers in Physics 9, 268 (1995)

/*!
 *  \brief Inverse hyperbolic cosine
 *  \param p_arg - input argument
 *  \return output value T type
 */
template<class T>
inline T acosh(T p_arg)
{
    return std::log(p_arg + std::sqrt(p_arg * p_arg - 1));
}

/*!
 *  \brief jacobian elliptic arctan function
 *  \param p_eival - the value of the elliptic integral \f$ u = \int_{0}^{v}(1 - m \sin^2 t)^{-1/2}dt \f$
 *  \param p_eimod - the modulus of the integral - m
 *  \return jacobian elliptic arctan function
 */
template<class T>
inline T arcsc(T p_eival, T p_eimod)
{
    std::size_t i, L;
    T A, B, BT, Y;

    A = 1.;
    B = p_eimod;
    Y = 1.0 / p_eival;
    L = 0;

    for(i = 0; i < MAX_TERMS; ++i)
    {
        BT = A * B;
        A = A + B;
        B = 2 * std::sqrt(BT);
        Y = Y - BT / Y;
        if(Y == 0){
            Y = std::sqrt(BT) * ERR_SMALL;}
        if(std::fabs(A - B) < (A * ERR_SMALL)){break;}
        L = 2 * L;
        if(Y < 0.){L++;}
    }
    if(Y < 0.){L++;}

    return ((std::atan(A / Y) + ADF_PI * L) / A);
}

/*!
 *  \brief Inverse hyperbolic sine
 *  \param p_arg - input argument
 *  \return output value T type
 */
template<class T>
inline T asinh(T p_arg)
{
    return std::log(p_arg + std::sqrt(p_arg * p_arg + 1));
}

/*!
 *  \brief This routine calculate jacobian elliptic functions using
 *         arithmetic-geometric mean(agm) method.
 *         Starting with \f$ a_0 = 1; b_0 = (1-m)^{1/2}; c_0 = m^{1/2} \f$
 *         we determine \f$ a_1, b_1, c_1; \ldots; a_n, b_n, c_n \f$
 *         with \f$ a_k = (a_{k-1} + b_{k-1})/2;  b_k = (a_{k-1} \times b_{k-1})^{1/2}; c_k = (a_{k-1} - b_{k-1})/2 \f$
 *         and we stop when \f$ c_n \approx 0 (in fact, when a_{n-1} = a_n) \f$
 *         This program calculates only the 3 copolar trio: sn(u|m); cn(u|m); dn(u|m).
 *         When the agm process stops, we calculates \f$ v_n = 2^{n} \times a_{n} \times u \times 180/\pi \f$ (in degrees)
 *         and then, \f$ v_{n-1}, \ldots, v_0 \f$ from the recurrence relation:
 *         \f$ \sin(2v_{k-1}-v_{k}) = (c_{k}/a_{k})\sin v_k \f$
 *         then,
 *         \f$ sn(u|m) = \sin v_0; cn(u|m) = \cos v_0 \f$ and \f$ dn(u|m) = (1 - m\sin^2 v_0)^{1/2} \f$
 *  \param p_eival - the value of the elliptic integral \f$ u = \int_{0}^{v}(1 - m \sin^2 t)^{-1/2}dt \f$
 *  \param p_eimod - the modulus of the integral - m
 *  \param sn - jacobian elliptic sine functions
 *  \param cn - jacobian elliptic cosine functions
 *  \param dn - jacobian elliptic differenc functions
 */
template<class T>
inline void ellip_funcs(T p_eival, T p_eimod, T &sn, T &cn, T &dn)
{
    int i, imax;
    std::array<T, MAX_TERMS> A, B, C;
    std::array<T, MAX_TERMS> P;

    p_eimod = p_eimod * p_eimod;

    A[0] = 1;
    B[0] = std::sqrt(1 - p_eimod);
    C[0] = std::sqrt(p_eimod);

    for(i = 1; i < MAX_TERMS; ++i)
    {
        A[i] = (A[i - 1] + B[i - 1]) / 2;
        B[i] = std::sqrt(A[i - 1] * B[i - 1]);
        C[i] = (A[i - 1] - B[i - 1]) / 2;

        if(C[i] < ERR_SMALL) { break; }
    }

    if(i == MAX_TERMS) {imax = i - 1;}
    else{imax = i;}

    P[imax] = std::pow(2, imax) * A[imax] * p_eival;
    for(i = imax; i > 0; --i)
    {
        P[i - 1] = (std::asin(C[i] * std::sin(P[i]) / A[i]) + P[i]) / 2;
    }

    sn = std::sin(P[0]);
    cn = std::cos(P[0]);
    dn = std::sqrt(1 - p_eimod * (sn * sn));
}

/*!
 *  \brief Algorithm for calc complete elliptic integral of the first kinds using AGM process.
 *         Start with three initialization values \f$(a_0, b_0, c_0), where
 *         \f$(a_0 - 1, b_0 - \sqrt{1-k^2}, c_0 - k)
 *         and on the i-th iteration cycle, the value will be as follows:
 *         \f$(
 *              a_i = (a_{i-1} + b_{i-1})/2,
 *              b_i = \sqrt{a_{i-1} + b_{i-1}},
 *              c_i = (a_{i-1} - b_{i-1})/2
 *            )
 *         The cycle ends when \f$( a_N = b_N ) within some tolerance \f$( \varepsilon ).
 *         And then the value of the integral will be equal to:
 *         \f$(
 *              E=\frac{\pi}{2a_N}
 *            )
 *          See Satinder S. Sidhu. Elliptic integrals and functions.
 *          https://aip.scitation.org/doi/pdf/10.1063/1.168529
 *  \param p_eimod - equal to k - the elliptic modulus
 *  \return value of the complete elliptic integral of the first kinds
 */
template<class T>
T ellip_integral(T p_eimod)
{
    int i;
    std::array<T, MAX_TERMS> A, B, C;

    A[0] = 1;
    B[0] = std::sqrt(1 - (p_eimod * p_eimod));
    C[0] = p_eimod;

    for(i = 1; i < MAX_TERMS; ++i)
    {
        A[i] = (A[i - 1] + B[i - 1]) / 2;
        B[i] = std::sqrt(A[i - 1] * B[i - 1]);
        C[i] = (A[i - 1] - B[i - 1]) / 2;

        if(C[i] < ERR_SMALL) {break;}
    }
    return ADF_PI/(2*A[i]);
}

/*!
 *  \brief Algorithm for calculate I_0 is the zeroth-order modified Bessel function of the first kind.
 *         Using in window functions finite impulse response filter design and spectral analysis.
 *         In general, when using the Kaiser-Bessel window function, the modified Bessel function is defined as:
 *         \f$(
 *              I_0(x) = \sum_{k=1}^{\infty} \bigg[ \frac{(\frac{x}{2})^k}{k!} \bigg]^{2}
 *            )
 *         See Rabiner L.R.-Theory and Application of Digital Signal Processing
 *  \param p_arg - the argument of Bessel function
 *  \return value of the modified Bessel function of the first kind
 */
template<class T>
T bessel_func_mod(T p_arg)
{
    std::size_t i, converge;
    T Iold, Inew, J, K;

    Iold = 1.;
    J = 1.;
    K = p_arg / 2.;
    converge = 0;

    for(i = 1; i < MAX_TERMS; ++i)
    {
        J *= K / i;
        Inew = Iold + (J * J);
        if((Inew - Iold) < ERR_SMALL){
            converge = 1;
            break;}
        Iold = Inew;
    }
    if(!converge){return 0.0;}
    return Inew;
}

}
#endif // ADVMATH_H
