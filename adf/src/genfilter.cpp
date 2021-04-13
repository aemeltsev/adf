#include "genfilter.hpp"
namespace adf {

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam, FilterSelect &sfilter, ApproxSelect &sapprox) noexcept
    :m_sfilter(sfilter)
    ,m_sapprox(sapprox)
{
    m_fparam = std::make_unique<FiltParam<T>>();
    m_fparam->gain_passband = std::move(fparam->gain_passband);
    m_fparam->gain_stopband = std::move(fparam->gain_stopband);
    m_fparam->freq_passband = std::move(fparam->freq_passband);
    m_fparam->freq_stopband = std::move(fparam->freq_stopband);
    m_fparam->fsamp = std::move(fparam->fsamp);
    m_fparam->gain = std::move(fparam->gain);
    //m_order = std::move(fparam->order);

}

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs() noexcept
{
    m_fparam = std::make_unique<FiltParam<T>>();
}

template<typename T>
void CalcFilterCoefs<T>::setFiltParam(std::pair<T, T>& g_passband,
        std::pair<T, T>& g_stopband,
        std::pair<T, T>& f_passband,
        std::pair<T, T>& f_stopband,
        T fsamp,
        T gain)
{
    m_fparam->gain_passband = std::move(g_passband);
    m_fparam->gain_stopband = std::move(g_stopband);
    m_fparam->freq_passband = std::move(f_passband);
    m_fparam->freq_stopband = std::move(f_stopband);
    m_fparam->fsamp = std::move(fsamp);
    m_fparam->gain = std::move(gain);
    /*m_fparam->order = std::move(order)*/;
}

template<typename T>
void CalcFilterCoefs<T>::setTypeFilter(FilterSelect& sfilter)
{
    m_sfilter = sfilter;
}

template<typename T>
void CalcFilterCoefs<T>::setApproxFilter(ApproxSelect& sapprox)
{
    m_sapprox = sapprox;
}

/**
 * @brief This the common value for all filter approximation methods
 *        \f$(\varepsilon_s / \varepsilon_p)
 *        \f$(\varepsilon_s = 10.0^{-0.1*a_s}-1) stopband gain adjustment factor and
 *        passband gain ratio \f$(\varepsilon_p = 10.0^{-0.1*a_p}-1)
 * @return Ratio value of the suppression \f$(R_s)dB/\f$(R_p)dB
 */
template<typename T>
T CalcFilterCoefs<T>::CommonKernel()
{
    return ((std::pow(10.0,-0.1*m_fparam->gain_stopband.first)-1)/
            (std::pow(10.0,-0.1*m_fparam->gain_passband.first)-1));
}

/**
 * @brief
 */
template<typename T>
void CalcFilterCoefs<T>::FilterOrder()
{
    T kernel, ratio,      /**< Internal consts */
      order,              /**< Order value */
      wp1, wp2, ws1, ws2; /**< Edge frequency variables */

    T ratio_const, kernel_const,          /**< TODO */
      ei_ratio_const, eiq_ratio_const,    /**< TODO */
      ei_kernel_const, eiq_kernel_const;

    wp1 = m_fparam->f_passband.first;
    wp2 = m_fparam->f_passband.second;
    ws1 = m_fparam->f_stopband.first;
    ws2 = m_fparam->f_stopband.second;

    switch(m_sfilter)
    {
    case FilterSelect::LPF:
        ratio = ws1/wp1;
        break;
    case FilterSelect::HPF:
        ratio = wp1/ws1;
        break;
    case FilterSelect::PBF:
        if(ws1 > (wp1 * wp2) / ws2)
        {
            ws2 = (wp1 * wp2) / ws1;
            m_fparam->f_stopband.first = ws2;
        }
        else
        {
            ws1 = (wp1 * wp2) / ws2;
            m_fparam->f_stopband.second = ws1;
        }
        ratio = (ws2 - ws1) / (wp2 - wp1);
        break;
    case FilterSelect::SBF:
        if(wp1 > (ws1 * ws2) / wp2)
        {
            wp2 = (ws1 * ws2) / wp1;
            m_fparam->f_passband.first = wp2;
        }
        else
        {
            wp1 = (ws1 * ws2) / wp2;
            m_fparam->f_passband.second = wp1;
        }
        ratio = (wp2 - wp1) / (ws2 - ws1);
        break;
    default: return /* error code */;
    }

    kernel = CommonKernel();

    switch(m_sapprox)
    {
    case ApproxSelect::BUTTER:
        order = std::log10(kernel)/(2 * std::log10(ratio));
        break;
    case ApproxSelect::CHEBY:
    case ApproxSelect::ICHEBY:
        order = acosh(std::sqrt(kernel))/acosh(ratio);
        break;
    case ApproxSelect::ELLIPT:
        ratio_const = 1/ratio;
        if(ratio_const > .9999) return /* error code */;
        kernel_const = 1/std::sqrt(kernel);
        if(ratio_const < 2e-8) return /* error code */;
        ei_ratio_const = ellip_integral(ratio_const);
        eiq_ratio_const = ellip_integral(std::sqrt(1-(ratio_const*ratio_const)));
        ei_kernel_const = ellip_integral(kernel_const);
        eiq_kernel_const = ellip_integral(std::sqrt(1-(kernel_const*kernel_const)));
        order = (ei_ratio_const * eiq_kernel_const)/(eiq_ratio_const * ei_kernel_const);
        break;
    default: return /* error code */;
    }
    if(order > 200) return /* error code */;
    m_order = static_cast<int16_t>(order);
}

/**
 * @brief 1. Calculate \f$( \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} )
 *        2. Calculate radius: \f$( R = \varepsilon^{-1/n} )
 *        3. In for cycle calculate stable function left-half-plane poles used:
 *           \f$(
 *                s_k=\omega_c\bigg[ -\sin\frac{(2K + 1)\pi}{2n} + j\cos\frac{(2K + 1)\pi}{2n} \bigg],~~~K=0,1,\cdots ,n-1
 *              )
 *            and angle:
 *            \f$(
 *                 \phi = \frac{\pi}{n} \times \frac{(2k+n+1)}{2}
 *               )
 */
template<typename T>
void CalcFilterCoefs<T>::ButterApprox()
{
    T epsilon, radius, /**< Ripple factor, radius of the circle*/
      theta, sigma, omega; /**< Real and image position in s-domain value \f$( s = \sigma + j\omega) */
    if(m_order <= 0) return ADF_Error(BadValue, "Error: Using bad value");

    epsilon = std::sqrt(std::pow(10.0, -0,1*m_fparam->g_passband.first) - 1);
    radius = std::pow(epsilon, -1.0/m_order);

    m_fparam->gain = 1.;
    int32_t a=0, b=0; /**< Counters */
    /**< Work with the odd order */
    if(m_order % 2){
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = radius;
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = radius;
    }
    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        /**< Calculate angle first, then real and imag pos */
        theta = ADF_PI*(2*m + m_order + 1) / (2*m_order);
        sigma = radius * std::cos(theta);
        omega = radius * std::sin(theta);

        /**< Set the quadratic coefs */
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sigma*sigma + omega*omega;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = -2*sigma;
        n_bcoefs[b++] = sigma*sigma + omega*omega;
    }
}

/**
 * @brief 1. Calculate \f$( \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} )
 *        2. Calculate radius(d): \f$( D = \frac{\sinh^{-1}(\varepsilon^{-1})}{n} )
 *        3. Calculate angle:
 *           \f$(
 *                \phi = \frac{\pi}{n} \times \frac{(2k+1)}{2}
 *              )
 *            And poles in the left half of the plane are given by
 *            \f$(
 *                 p_k = -\sinh(D)\sin(\phi_k)+j\cosh(D)\cos(\phi_k)~~~k = 1, 2, 3, \cdots, n
 *               )
 */
template<typename T>
void CalcFilterCoefs<T>::ChebyApprox()
{
    T epsilon, d, /**< Ripple factor, \sigma axis radius */
      phi, sigma, omega; /**< Angle, real and image position in s-domain value \f$( s = \sigma + j\omega) */

    if(m_order <= 0) return ADF_Error(BadValue, "Error: Using bad value");

    epsilon = std::sqrt(std::pow(10.0, -0,1*m_fparam->g_passband.first) - 1);
    d = asinh(1/epsilon) / m_order;

    int32_t a=0, b=0; /**< Counters */
    /**< Work with the odd order */
    if(m_order % 2){
        m_fparam->gain = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sinh(d);
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = sinh(d);
    }
    else{
        m_fparam->gain = std::pow(10.,0.05 * m_fparam->g_passband.first);
    }

    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        /**< Calculate angle first, then real and imag pos */
        phi = ADF_PI*(2*m + 1) / (2*m_order);
        sigma = -1 * sinh(d) * std::sin(phi);
        omega = cosh(d) * std::cos(phi);

        /**< Set the quadratic coefs */
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sigma*sigma + omega*omega;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = -2*sigma;
        n_bcoefs[b++] = sigma*sigma + omega*omega;
    }
}

/**
 * @brief 1. Calculate \f$( \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} )
 *        TODO
 */
template<typename T>
void CalcFilterCoefs<T>::ElliptApprox()
{
    T epsilon,                          /**< Ripple factor */
      ratio,                            /**< Check type filter frequency value */
      sp, cp, dp,                       /**< Sn cn dn Jacobi elliptic functions */
      sn, cn, dn,
      kernel,                           /**< Internal const */
      ratio_const, kernel_const,        /**< Internal const */
      ei_ratio_const, ei_kernel_const,  /**< Complete elliptic integral from kernel and ratio const */
      odd,                              /**< Check order value */
      vo, fm,                           /**< Internal const */
      sigma, omega, zero;               /**< Pole and zero */

    if(m_order <= 0) return ADF_Error(BadValue, "Error: Using bad value");

    epsilon = std::sqrt(std::pow(10., -0,1*m_fparam->g_passband.first) - 1);

    switch (m_sfilter)
    {
    case FilterSelect::LPF:
        ratio = m_fparam->f_stopband.first / m_fparam->f_passband.first;
        break;
    case FilterSelect::HPF:
        ratio = m_fparam->f_passband.first / m_fparam->f_stopband.first;
        break;
    case FilterSelect::PBF:
        ratio = (m_fparam->f_stopband.second - m_fparam->f_stopband.first) /
                (m_fparam->f_passband.second - m_fparam->f_passband.first);
        break;
    case FilterSelect::SBF:
        ratio = (m_fparam->f_passband.second - m_fparam->f_passband.first) /
                (m_fparam->f_stopband.second - m_fparam->f_stopband.first);
        break;
    default: return /* error code */;
    }

    kernel = CommonKernel();
    ratio_const = 1/ratio;
    kernel_const = 1/std::sqrt(kernel);
    ei_ratio_const = ellip_integral(ratio_const);
    ei_kernel_const = ellip_integral(kernel_const);
    vo = (ei_ratio_const / (ei_kernel_const * m_order)) * arcsc(1/epsilon, kernel_const);
    ellip_funcs(vo, std::sqrt(1-(ratio_const*ratio_const)), sp, cp, dp);

    int32_t a=0, b=0; /**< Counters */
    odd = m_order % 2;
    if(odd){
        m_fparam->gain = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sp * cp /(1 - sp*sp);
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = sp * cp /(1 - sp*sp);
    }
    else{
        m_fparam->gain = std::pow(10., 0.05 * m_fparam->g_passband.first);
    }

    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        /**< Calc angle first, then real and imag pos */
        fm = ei_ratio_const * (2*m + 1 + odd) / m_order;
        ellip_funcs(fm, ratio_const, sn, cn, dn);

        /**< Calculated real and imag coordinates of poles */
        sigma = -1 * cn*dn*sp*cp /(1 - dn*dn*sp*sp);
        omega = sn*dp / (1 - dn*dn*sp*sp);

        /**< Calculated the zero location */
        zero = 1 / (ratio_const * sn);

        /**< Set the quadratic coefs */
        n_acoefs[a++] = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = zero*zero;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = -2*sigma;
        n_bcoefs[b++] = sigma*sigma + omega*omega;

        /**< Update the gain */
        m_fparam->gain *=
                ((sigma*sigma + omega*omega)/(zero*zero));
    }
}

/**
 * @btief TODO
 */
template<typename T>
void CalcFilterCoefs<T>::IChebyApprox()
{
    T epsilon, d,         /**< Ripple factor */
      mag_inv,            /**< Using for inverse magnitude value */
      phi, sigma, omega,  /**< Real and image position in s-domain value \f$( s = \sigma + j\omega) */
      zero;               /**< Zero value */

    if(m_order <= 0) return ADF_Error(BadValue, "Error: Using bad value");

    epsilon = std::sqrt(std::pow(10.0, -0,1*m_fparam->g_stopband.first) - 1.);
    d = asinh(1/epsilon) / m_order;

    m_fparam->gain = 1.;
    int32_t a=0, b=0; /**< Counters */

    /**< work with the odd order */
    if(m_order % 2){
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 1 / sinh(d);
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = 1 / sinh(d);
    }

    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        /**< Calculate angle first, then real and imag position */
        phi = ADF_PI*(2*m + 1) / (2*m_order);
        sigma = -1 * sinh(d) * std::sin(phi);
        omega = cosh(d) * std::cos(phi);

        /**< Calculate inverse of pole location */
        mag_inv = omega*omega + sigma*sigma;
        omega = -1 * omega/mag_inv;
        sigma = sigma/mag_inv;

        /**< Calculate the zero location */
        zero = 1. / std::cos(phi);

        /**< Set the quadratic coefs */
        mag_inv = omega*omega + sigma*sigma;
        n_acoefs[a++] = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = zero * zero;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = -2*sigma;
        n_bcoefs[b++] = mag_inv;

        /**< Update the gain */
        m_fparam->gain *= (mag_inv / (zero * zero));
    }
}

template<typename T>
void CalcFilterCoefs<T>::NormalCoefs()
{
    if(m_order <= 0) return ADF_Error(BadValue, "Error: Using bad value");
    int32_t number_coefs = 3 * (m_order + 1) / 2;
    n_acoefs.reserve(number_coefs);
    n_bcoefs.reserve(number_coefs);

    switch (m_sapprox)
    {
    case ApproxSelect::BUTTER:
        ButterApprox();
        break;
    case ApproxSelect::CHEBY:
        ChebyApprox();
        break;
    case ApproxSelect::ELLIPT:
        ElliptApprox();
        break;
    case ApproxSelect::ICHEBY:
        IChebyApprox();
        break;
    default: return ADF_Error(BadFilter, "Error: Bad approximation value");
    }

}

template<typename T>
void CalcFilterCoefs<T>::BSCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::BPCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::HPCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::LPCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::UnnormCoefs()
{

}

}
