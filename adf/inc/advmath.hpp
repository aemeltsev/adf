/**
  * License text...
*/
#ifndef ADVMATH_H
#define ADVMATH_H
#include <cstdint>
#include <cmath>
#include <limits>
#include <cfloat>
#include <vector>
#include <array>

#include "common.hpp"

namespace adf {

/**
 *  @brief
 *  @param
 *  @return
 */
template<class T>
inline T acosh(T p_arg)
{
    return std::log(p_arg + std::sqrt(p_arg*p_arg-1));
}

/**
 *  @brief
 *  @param
 *  @param
 *  @param
 *  @return
 */
template<class T>
inline T arcsc(T p_eival, T p_eimod, int32_t p_prec)
{
    int i,L;
    T A, B, BT, Y;

    A = 1.;
    B = p_eimod;
    Y = 1.0/p_eival;
    L = 0;

    for(i=0; i<p_prec; i++)
    {
        BT = A * B;
        A = A + B;
        B = 2 * std::sqrt(BT);
        Y = Y - BT / Y;
        if(Y == 0){
            Y = std::sqrt(BT)*FLT_MIN;}
        if(std::fabs(A-B)<(A*FLT_MIN)){break;}
        L = 2 * L;vector>
        if(Y < 0.){L++;}
    }
    if(Y < 0.){L++;}

    return ((std::atan(A/Y) + M_PI * L) / A);
}

/**
 *  @brief
 *  @param
 *  @return
 */
template<class T>
inline T asinh(T p_arg)
{
    return std::log(p_arg + std::sqrt(p_arg*p_arg+1));
}

/**
 *  @brief
 *  @param
 *  @param
 *  @param
 *  @param
 *  @param
 *  @param
 *  @return
 */
template<class T>
inline void ellip_funcs(T p_eival, T p_eimod, T &sn, T &cn, T &dn, int32_t p_prec)
{
    int16_t i, imax;
    std::vector<T> A, B, C;
    std::vector<T> P;

    p_eimod=p_eimod*p_eimod;

    A[0]=1;
    B[0]=std::sqrt(1-p_eimod);
    C[0]=std::sqrt(p_eimod);

    for(i=1; i<p_prec; i++)
    {
        A[i]=(A[i-1]+B[i-1])/2;
        B[i]=std::sqrt(A[i-1]*B[i-1]);
        C[i]=(A[i-1]-B[i-1])/2;

        if(C[i]<FLT_MIN){
            break;}
    }

    if(i == p_prec){
        imax=i-1;}
    else{
        imax=i;}

    P[imax]=std::pow(2,imax)*A[imax]*p_eival;
    for(i=imax; i>0; i--)
    {
        P[i-1]=(std::asin(C[i]*std::sin(P[i])/A[i])+P[i])/2;
    }

    sn=std::sin(P[0]);
    cn=std::cos(P[0]);
    dn=std::sqrt(1-p_eimod*(sn*sn));
}

/**
 *  @brief Algorithm for calc complete elliptic integral of the first kinds using AGM process.
 *         Start with three initialization values \f$(a_0, b_0, c_0), where
 *         \f$(a_0 - 1, b_0 - \sqrt{1-k^2}, c_0 - k)
 *         and on the i-th iteration cycle, the value will be as follows:
 *         \f$(
 *              a_i = (a_{i-1} + b_{i-1})/2,
 *              b_i = \sqrt{a_{i-1} + b_{i-1}},
 *              c_i = (a_{i-1} - b_{i-1})/2
 *            )
 *         Еhe cycle ends when \f$( a_N = b_N ) within some tolerance \f$( \varepsilon ).
 *         And then the value of the integral will be equal to:
 *         \f$(
 *              E=\frac{\pi}{2a_N}
 *            )
 *  @param p_eimod - equal to k - the elliptic modulus
 *  @return value of the complete elliptic integral of the first kinds
 */
template<class T>
T ellip_integral(T p_eimod)
{
    int16_t i;
    std::array<T, MAX_TERMS> A, B, C;

    A[0]=1;
    B[0]=std::sqrt(1-(p_eimod*p_eimod));
    C[0]=p_eimod;

    for(i=1; i<MAX_TERMS; i++)
    {
        A[i]=(A[i-1]+B[i-1])/2;
        B[i]=std::sqrt(A[i-1]*B[i-1]);
        C[i]=(A[i-1]-B[i-1])/2;

        if(C[i]<ERR_SMALL){
            break;}
    }
    return ADF_PI/(2*A[i]);
}

/**
 *  @brief
 *  @param
 *  @param
 *  @return
 */
template<class T>
T bessel_func_mod(T p_arg, int32_t p_prec)
{
    int32_t i, converge;
    T Iold, Inew, J, K;

    Iold=1.;
    J=1.;
    K=p_arg/2.;
    converge=0;

    for(i=1; i<p_prec; i++)
    {
        J *= K/i;
        Inew = Iold+(J*J);
        if((Inew-Iold)<FLT_MIN){
            converge=1;
            break;}
        Iold=Inew;
    }
    if(!converge){return 0.0;}
    return Inew;
}

}

#endif // ADVMATH_H
