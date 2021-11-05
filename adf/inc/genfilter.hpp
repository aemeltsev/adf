#ifndef GENFILTER_H
#define GENFILTER_H
#include <cstdint>
#include <memory>
#include <cmath>
#include <vector>

#include "advmath.hpp"
#include "complex.hpp"
#include "base.hpp"

namespace adf {

/**
 * @brief The FilterType enum - the enumeration class for the type filter select
 */
enum class FilterType
{
    LPF=1, //Low-pass filter
    HPF,   //High-pass filter
    PBF,   //Band-pass filter
    SBF    //Band-stop filter
};

/**
 * @brief The ApproxType enum - the enumeration class for select the approximation methods
 */
enum class ApproxType
{
    BUTTER=1, //Butterworth approximation
    CHEBY,    //Chebyshev approximation
    ICHEBY,   //Inverse Chebyshev approximation
    ELLIPT    //elliptic approximation
};

/**
 * @brief Input specification for filter parameters,
 *        these include the passband and stopband edge frequency
 *        (in freq_passband first value - frequency lower, second value - frequency upper,
 *        in freq_stopband the first value - frequency lower, the second value - frequency upper)
 *        and gains(in gain_passband and gain_stopband values)
 *        used double type by default
 * @param gain_passband
 * @param gain_stopband
 * @param freq_passband
 * @param freq_stopband
 * @param fsamp - sampe frequency for digital filter
 * @param gain - gain multiplier
 * @param order - order, length of filter
 */
template<typename T=double>
struct FiltParam
{
  std::pair<T, T> gain_passband;
  std::pair<T, T> gain_stopband;
  std::pair<T, T> freq_passband;
  std::pair<T, T> freq_stopband;
  T fsamp;
  T gain;
};

/**
 * @brief
 */
template<typename T=double>
class CalcFilterCoefs
{
private:
    FiltParam<T> m_fparam;
    FilterType m_sfilter;
    ApproxType m_sapprox;
    std::size_t m_order = 0;
    std::size_t m_gain = 0;
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */

protected:
    T CommonKernel();
    T FreqNorm();
    void FilterOrder();
    bool FillZeroCoeffs(std::vector<T>& avec, std::vector<T>& bvec, std::size_t num);

    void NormalCoefs();
    void ButterApprox();
    void ChebyApprox();
    void ElliptApprox();
    void IChebyApprox();

    void UnnormCoefs();
    void BSCoefsUnnorm(T un_bandwith, T un_centrfreq);
    void BPCoefsUnnorm(T un_bandwith, T un_centrfreq);
    void HPCoefsUnnorm(T freq);
    void LPCoefsUnnorm(T freq);

public:
    explicit CalcFilterCoefs(const FiltParam<T> &fparam, const FilterType &fselect, const ApproxType &sapprox) noexcept
        :m_sfilter(fselect)
        ,m_sapprox(sapprox)
    {
        m_fparam.gain_passband = std::move(fparam.gain_passband);
        m_fparam.gain_stopband = std::move(fparam.gain_stopband);
        m_fparam.freq_passband = std::move(fparam.freq_passband);
        m_fparam.freq_stopband = std::move(fparam.freq_stopband);
        m_fparam.fsamp = std::move(fparam.fsamp);
        m_fparam.gain = std::move(fparam.gain);
    }

    CalcFilterCoefs() noexcept
    {
        m_order = 1;
        m_gain = 1;
    }

    /**
     * @brief setFiltParam Filling data fields of the
     * @param g_passband - the passband gain ripple
     * @param g_stopband - the stopband gain ripple
     * @param f_passband - the passband edge frequency
     * @param f_stopband - the stopband edge frequency
     * @param fsamp - sample frequency
     * @param gain - gain miltiplier
     */
    void setFiltParam(
            const std::pair<T, T> &g_passband,
            const std::pair<T, T> &g_stopband,
            const std::pair<T, T> &f_passband,
            const std::pair<T, T> &f_stopband,
            const T fsamp,
            const T gain)
    {

        m_fparam.gain_passband = std::move(g_passband);
        m_fparam.gain_stopband = std::move(g_stopband);
        m_fparam.freq_passband = std::move(f_passband);
        m_fparam.freq_stopband = std::move(f_stopband);
        m_fparam.fsamp = fsamp;
        m_fparam.gain = gain;
    }

    void setTypeFilter(const FilterType& sfilter)
    {
        m_sfilter = sfilter;
    }

    void setApproxFilter(const ApproxType& sapprox)
    {
        m_sapprox = sapprox;
    }

    // Only for test using ----------------------------------------------------------
    void makeFilterOrder(){FilterOrder();}
    void makeNormalCoefs(){NormalCoefs();}
    void makeUnnormalCoefs(){UnnormCoefs();}

    std::size_t getFilterOrder() const
    {
        return m_order;
    }

    ApproxType getApproxType()
    {
        return m_sapprox;
    }

    FilterType getFilterType()
    {
        return m_sfilter;
    }

    std::vector<T> normACoefs()
    {
        return n_acoefs;
    }

    std::size_t normACoefsSize()
    {
        return n_acoefs.size();
    }

    std::vector<T> normBCoefs()
    {
        return n_bcoefs;
    }

    std::size_t normBCoefsSize()
    {
        return n_bcoefs.size();
    }

    std::vector<T> unnormACoefs()
    {
        return un_acoefs;
    }

    std::size_t unnormACoefsSize()
    {
        return un_acoefs.size();
    }

    std::vector<T> unnormBCoefs()
    {
        return un_bcoefs;
    }

    std::size_t unnormBCoefsSize()
    {
        return un_bcoefs.size();
    }
};

/**
 * @brief This the common value for all filter approximation methods
 *        \f$\varepsilon_s / \varepsilon_p\f$
 *        \f$\varepsilon_s = 10.0^{-0.1*a_s}-1\f$ stopband gain adjustment factor and
 *        passband gain ratio \f$\varepsilon_p = 10.0^{-0.1*a_p}-1\f$
 * @return Ratio value of the suppression \f$(R_s)dB/\f$(R_p)dB\f$
 */
template<typename T>
T CalcFilterCoefs<T>::CommonKernel()
{
    if(m_fparam.gain_stopband.first <=ADF_GAIN_STOP || m_fparam.gain_passband.first >= ADF_GAIN_PASS)
    {
        throw std::invalid_argument(ADF_ERROR("Zero or negative gain value"));
    }
    return ((std::pow(10.0,-0.1*m_fparam.gain_stopband.first)-1)/
            (std::pow(10.0,-0.1*m_fparam.gain_passband.first)-1));
}

/**
 * @brief The normalization to relative of the passband cutoff frequency
 *        \f$\omega = \fract{f_1}{f_2}\f$
 *        The part of filter order calculation
 * @return Ratio value of the normalization by frequency
 */
template <typename T>
T CalcFilterCoefs<T>::FreqNorm()
{
    /**< Return ratio value */
    T ratio;

    /**< Edge frequency variables */
    auto&& wp1 = m_fparam.freq_passband.first; /**<pb freq lower */
    auto&& wp2 = m_fparam.freq_passband.second; /**<pb freq upper */
    auto&& ws1 = m_fparam.freq_stopband.first; /**<sb freq lower */
    auto&& ws2 = m_fparam.freq_stopband.second; /**<sb freq upper */

    switch(m_sfilter)
    {
    case FilterType::LPF:
        ratio = ws1/wp1;
        break;
    case FilterType::HPF:
        ratio = wp1/ws1;
        break;
    case FilterType::PBF:
        if(ws1 > (wp1 * wp2) / ws2)
        {
            ws2 = (wp1 * wp2) / ws1;
            m_fparam.freq_stopband.second = ws2;
        }
        else
        {
            ws1 = (wp1 * wp2) / ws2;
            m_fparam.freq_stopband.first = ws1;
        }
        ratio = (ws2 - ws1) / (wp2 - wp1);
        break;
    case FilterType::SBF:
        if(wp1 > (ws1 * ws2) / wp2)
        {
            wp2 = (ws1 * ws2) / wp1;
            m_fparam.freq_passband.second = wp2;
        }
        else
        {
            wp1 = (ws1 * ws2) / wp2;
            m_fparam.freq_passband.first = wp1;
        }
        ratio = (wp2 - wp1) / (ws2 - ws1);
        break;
    default:
        ratio=0.;
        //throw std::invalid_argument(ADF_ERROR("Use undefine filter type"));
        break;
    }

    return ratio;
}

/**
 * @brief The order of the polynomial for the approximation,
 *        as defined in the filter specification, which is the order of the filter.
 *        As example Butterworth:
 *        \f[ n \geq \frac{\log(\frac{\varepsilon_s^2}{\varepsilon_p^2})}{2\log(\frac{\omega_s}{\omega_p})} ]\f
 * @return Order of the polynom value
 */
template<typename T>
void CalcFilterCoefs<T>::FilterOrder()
{
    T ratio_const, kernel_const,          /**< Temp values, */
      ei_ratio_const, eiq_ratio_const,    /**< for elliptic approximation */
      ei_kernel_const, eiq_kernel_const,
      order=0.;

    auto kernel = CommonKernel();
    auto ratio = FreqNorm();

    switch(m_sapprox)
    {
    case ApproxType::BUTTER:
        order = std::log10(kernel)/(2 * std::log10(ratio));
        break;
    case ApproxType::CHEBY:
    case ApproxType::ICHEBY:
        order = acosh(std::sqrt(kernel))/acosh(ratio);
        break;
    case ApproxType::ELLIPT:
        ratio_const = 1/ratio;

        if((ratio_const > .9999)||(ratio_const < 2e-8))
        {
            m_order = static_cast<std::size_t>(order);
            throw std::range_error(ADF_ERROR("The value to out of range"));
        }

        kernel_const = 1/std::sqrt(kernel);

        ei_ratio_const = ellip_integral(ratio_const);
        eiq_ratio_const = ellip_integral(std::sqrt(1-(ratio_const*ratio_const)));
        ei_kernel_const = ellip_integral(kernel_const);
        eiq_kernel_const = ellip_integral(std::sqrt(1-(kernel_const*kernel_const)));
        order = (ei_ratio_const * eiq_kernel_const)/(eiq_ratio_const * ei_kernel_const);
        break;
    }

    if(order > 200.)
    {
        m_order = static_cast<std::size_t>(order);
        throw std::range_error(ADF_ERROR("The filter order very large size"));
    }
    m_order = static_cast<std::size_t>(std::ceil(order));
}

/**
 * @brief FillZeroCoeffs - Filling normalized vectors with default values
 * @param avec input reference to vector of "a" coefficients for filling
 * @param bvec input reference to vector of "b" coefficients for filling
 * @param num the order
 * @return bool variable, success if order value not out of range
 */
template<typename T>
bool CalcFilterCoefs<T>::FillZeroCoeffs(std::vector<T>& avec, std::vector<T>& bvec, std::size_t num)
{
    bool success = false;
    if(num > 0 && num < MAX_TERMS)
    {
        success = true;
        std::size_t number_coeffs = (3 * (num + 1)) / 2;
        avec.reserve(number_coeffs);
        bvec.reserve(number_coeffs);

        for(std::size_t ind=0; ind<number_coeffs; ++ind)
        {
            avec.push_back(0);
            bvec.push_back(0);
        }
    }
    return success;
};

/**
 * @brief NormalCoefs - Calculation and filling the vectors of coefficients depending on the approximation method
 */
template<typename T>
void CalcFilterCoefs<T>::NormalCoefs()
{
    if(FillZeroCoeffs(n_acoefs, n_bcoefs, m_order))
    {
        switch (m_sapprox)
        {
        case ApproxType::BUTTER:
            ButterApprox();
            break;
        case ApproxType::CHEBY:
            ChebyApprox();
            break;
        case ApproxType::ELLIPT:
            ElliptApprox();
            break;
        case ApproxType::ICHEBY:
            IChebyApprox();
            break;
        }
    }
    else
    {
        throw std::range_error(ADF_ERROR("The order value to out of range"));
    }
}

/**
 * @brief 1. Calculate \f$( \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} )\f$
 *        2. Calculate radius: \f$( R = \varepsilon^{-1/n} )\f$
 *        3. In for cycle calculate stable function left-half-plane poles used:
 *           \f$(
 *                s_k=\omega_c\bigg[ -\sin\frac{(2K + 1)\pi}{2n} + j\cos\frac{(2K + 1)\pi}{2n} \bigg],~~~K=0,1,\cdots ,n-1
 *              )\f$
 *            and angle:
 *            \f$(
 *                 \phi = \frac{\pi}{n} \times \frac{(2k+n+1)}{2}
 *               )\f$
 */
template<typename T>
void CalcFilterCoefs<T>::ButterApprox()
{
    //Determine ripple factor
    auto epsilon = std::sqrt(std::pow(10.0, -0.1*m_fparam.gain_passband.first) - 1);

    //Determine the Butterworth radius
    auto radius = std::pow(epsilon, -1.0/m_order);

    //Default gain
    m_fparam.gain = 1.;

    //Counters
    int32_t a=0, b=0;

    //Work with the odd order
    if(m_order % 2){
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = radius;
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = radius;
    }
    /**< Other all quadratic terms,  */
    for(std::size_t m=0; m<m_order/2; m++)
    {
        /**< First determine the angle,
         *  and then the position of the complex pole,
         *  its real and imaginary values. */
        auto theta = ADF_PI*(2*m + m_order + 1) / (2*m_order);
        auto sigma = radius * std::cos(theta);
        auto omega = radius * std::sin(theta);

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
 * @brief 1. Calculate \f$( \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} )\f$
 *        2. Calculate radius(d): \f$( D = \frac{\sinh^{-1}(\varepsilon^{-1})}{n} )\f$
 *        3. Calculate angle:
 *           \f$(
 *                \phi = \frac{\pi}{n} \times \frac{(2k+1)}{2}
 *              )\f$
 *            And poles in the left half of the plane are given by
 *            \f$(
 *                 p_k = -\sinh(D)\sin(\phi_k)+j\cosh(D)\cos(\phi_k)~~~k = 1, 2, 3, \cdots, n
 *               )\f$
 */
template<typename T>
void CalcFilterCoefs<T>::ChebyApprox()
{
    //Determine ripple factor
    auto epsilon = std::sqrt(std::pow(10.0, -0.1*m_fparam.gain_passband.first) - 1);

    //Determine minor axis radius of the ellipse
    auto d = asinh(1/epsilon) / m_order;

    // Counters
    int32_t a=0, b=0;

    /**< Work with the odd order */
    if(m_order % 2){
        m_fparam.gain = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sinh(d);
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = sinh(d);
    }
    else{
        m_fparam.gain = std::pow(10., 0.05*m_fparam.gain_passband.first);
    }

    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        /**< First determine the angle,
         *  and then the position of the complex pole,
         *  its real and imaginary values. */
        auto phi = ADF_PI*(2*m + 1) / (2*m_order);
        auto sigma = -1 * sinh(d) * std::sin(phi);
        auto omega = cosh(d) * std::cos(phi);

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
 * @brief 1. Calculate \f$ \varepsilon = \sqrt[]{\left( 10^{^{A_p}/_{10}}-1 \right)} \f$
 *        2. Check type filter for calculate the normalized cutoff frequency ratio
 *        3. Calculate kernel ratio, and temp variables
 *           \f$ v_0 = -\frac{j}{NK_1}sn^{-1} \bigg( \frac{j}{\varepsilon_p}, k_1 \bigg) \f$
 *        4. TODO
 */
template<typename T>
void CalcFilterCoefs<T>::ElliptApprox()
{
      T ratio,                            /**< Check type filter frequency value */
      sp, cp, dp,                       /**< Sn cn dn Jacobi elliptic functions */
      sn, cn, dn;

    // The attenuation unevenness ratio in passband - \f$ \varepsilon \f$ (Ripple factor)
    auto epsilon = std::sqrt(std::pow(10., -0.1*m_fparam.gain_passband.first) - 1);

    // Normalized cutoff fre\f$quency ratio
    switch (m_sfilter)
    {
    case FilterType::LPF:
        ratio = m_fparam.freq_stopband.first / m_fparam.freq_passband.first;
        break;
    case FilterType::HPF:
        ratio = m_fparam.freq_passband.first / m_fparam.freq_stopband.first;
        break;
    case FilterType::PBF:
        ratio = (m_fparam.freq_stopband.second - m_fparam.freq_stopband.first) /
                (m_fparam.freq_passband.second - m_fparam.freq_passband.first);
        break;
    case FilterType::SBF:
        ratio = (m_fparam.freq_passband.second - m_fparam.freq_passband.first) /
                (m_fparam.freq_stopband.second - m_fparam.freq_stopband.first);
        break;
    default:
        ratio=0.;
        throw std::invalid_argument(ADF_ERROR("Use undefine filter type"));
        break;
    }

    /* Dimopoulos H.G.-Analog Electronic Filters.Theory, Design and Synthesis.
     *  pp. 180. Specifications and the Order of the Elliptic Approximation */
    auto kernel = CommonKernel();
    auto ratio_const = 1/ratio;
    auto kernel_const = 1/std::sqrt(kernel);

    //The complete elliptic integrals of the modules ratio and kernel
    auto ei_ratio_const = ellip_integral(ratio_const); //K(k)
    auto ei_kernel_const = ellip_integral(kernel_const); //K(g)

    //Variable vo used in the calculation of the pole and zero locations
    auto vo = (ei_ratio_const / (ei_kernel_const * m_order)) * arcsc(1/epsilon, kernel_const); //N(d)
    ellip_funcs(vo, std::sqrt(1-(ratio_const*ratio_const)), sp, cp, dp);

    // Counters
    int32_t a=0, b=0;

    //Check odd filter order value
    auto odd = m_order % 2;
    if(odd){
        m_fparam.gain = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = sp * cp /(1 - sp*sp);
        n_bcoefs[b++] = 0.;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = sp * cp /(1 - sp*sp);
    }
    else{
        m_fparam.gain = std::pow(10., 0.05 * m_fparam.gain_passband.first);
    }

    /**< Other all quadratic terms */
    for(int32_t m=0; m<m_order/2; m++)
    {
        //Define fm variable using in calculation of the real and image parts
        auto fm = ei_ratio_const * (2*m + 1 + odd) / m_order;
        ellip_funcs(fm, ratio_const, sn, cn, dn);

        //Calculated real and imag coordinates of poles
        auto sigma = -1 * cn*dn*sp*cp /(1 - dn*dn*sp*sp);
        auto omega = sn*dp / (1 - dn*dn*sp*sp);

        //Calculated the zero location
        auto zero = 1 / (ratio_const * sn);

        /**< Set the quadratic coefs */
        n_acoefs[a++] = 1.;
        n_acoefs[a++] = 0.;
        n_acoefs[a++] = zero*zero;
        n_bcoefs[b++] = 1.;
        n_bcoefs[b++] = -2*sigma;
        n_bcoefs[b++] = sigma*sigma + omega*omega;

        /**< Update the gain */
        m_fparam.gain *=((sigma*sigma + omega*omega)/(zero*zero));
    }
}

/**
 * @brief IChebyApprox - Type II Chebyshev filter also inverse Chebyshev filter
 *        See details Dimopoulos H.G.-Analog Electronic Filters.Theory, Design and Synthesis.
 *        pp. 110
 *        The Inverse Chebyshev Approximation
 */
template<typename T>
void CalcFilterCoefs<T>::IChebyApprox()
{
    T  mag_inv,            /**< Using for inverse magnitude value */
      phi, sigma, omega,  /**< Real and image position in s-domain value \f$( s = \sigma + j\omega) */
      zero;               /**< Zero value */

    auto epsilon = std::sqrt(std::pow(10.0, -0.1*m_fparam.gain_stopband.first) - 1.);
    auto d = asinh(1/epsilon) / m_order;

    m_fparam.gain = 1.;

    //Counters
    int32_t a=0, b=0;

    //For the odd order - first order pole on the negative real axis */
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
        //Calculate angle
        auto phi = ADF_PI*(2*m + 1) / (2*m_order);

        //Next the pole location values must be inverted */
        auto isigma = -1 * sinh(d) * std::sin(phi);
        auto iomega = cosh(d) * std::cos(phi);

        //And calculate final pole locations
        mag_inv = iomega*iomega + isigma*isigma;
        omega = -1 * iomega/mag_inv;
        sigma = isigma/mag_inv;

        //Calculate the zero location
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
        m_fparam.gain *= (mag_inv / (zero * zero));
    }
}

/**
 * @brief
 */
template<typename T>
void CalcFilterCoefs<T>::UnnormCoefs()
{
    T freq, // Stored value for unnormaliztion frequency
      BW,   /* For transformation to bandstop or bandpass type */
      Wo;

    if(m_sapprox == ApproxType::BUTTER
       || m_sapprox == ApproxType::CHEBY
       || m_sapprox == ApproxType::ELLIPT)
    {
        freq = m_fparam.freq_passband.first;
        Wo = std::sqrt(m_fparam.freq_passband.first * m_fparam.freq_passband.second);
        BW = m_fparam.freq_passband.second - m_fparam.freq_passband.first;
    }
    else if(m_sapprox == ApproxType::ICHEBY)
    {
        freq = m_fparam.freq_stopband.first;
        Wo = std::sqrt(m_fparam.freq_stopband.first * m_fparam.freq_stopband.second);
        BW = m_fparam.freq_passband.second - m_fparam.freq_passband.first;
    }

    switch (m_sfilter)
    {
    case FilterType::LPF:
        LPCoefsUnnorm(freq);
        break;
    case FilterType::HPF:
        HPCoefsUnnorm(freq);
        break;
    case FilterType::PBF:
        BPCoefsUnnorm(BW, Wo);
        break;
    case FilterType::SBF:
        BSCoefsUnnorm(BW, Wo);
        break;
    }
}

/**
 * @brief
 */
template<typename T>
void CalcFilterCoefs<T>::BSCoefsUnnorm(T un_bandwith, T un_centrfreq)
{
    std::size_t origin_qd_count,                         /**< Original number of quads values */
                origin_order;                            /**< Original order */
    std::size_t origin_coef, new_coef, pos_start;        /**< Counters */
    std::size_t size_coef;                               /**< Size vector value */
    complex<T> A, B, C, D, E;                            /**< Temp complex value */

    /**
      Store the original number of the order,
      new order will be twice large to original
    */
    origin_order = m_order;
    origin_qd_count = (origin_order + 1) / 2;
    m_order = origin_order * 2;
    /**<  */
    size_coef = 3*origin_order;

    if(FillZeroCoeffs(un_acoefs, un_bcoefs, size_coef))
    {
        /**< If original order is odd convert first order factor to quadratic,
        *  pos_start indicate start point for loop */
        if(origin_order % 2)
        {
            m_gain *= (n_acoefs[2] / n_bcoefs[2]);
            un_acoefs[0] = 1.;
            un_acoefs[1] = un_bandwith * n_acoefs[1] / n_acoefs[2];
            un_acoefs[2] = un_centrfreq * un_centrfreq;
            un_bcoefs[0] = 1.;
            un_bcoefs[1] = un_bandwith * n_bcoefs[1] / n_bcoefs[2];
            un_bcoefs[3] = un_centrfreq * un_centrfreq;
            pos_start = 1;
        }
        else
        {
            pos_start = 0;
        }

        /**<  */
        for(std::size_t qd_count = pos_start; qd_count < origin_qd_count; qd_count++)
        {
            origin_coef = qd_count * 3;
            new_coef = qd_count * 6 - pos_start * 3;
            m_gain *= (n_acoefs[origin_coef+2] / n_bcoefs[origin_coef+2]);

            if(n_acoefs[origin_coef] == 0)
            {
                un_acoefs[new_coef] = 1.;
                un_acoefs[new_coef+1] = 0.;
                un_acoefs[new_coef+2] = un_centrfreq * un_centrfreq;
                un_acoefs[new_coef+3] = 1.;
                un_acoefs[new_coef+4] = 0.;
                un_acoefs[new_coef+5] = un_centrfreq * un_centrfreq;
            }
            /**<  */
            else
            {
                /**< Convert coefficients to complex, then factorization */
                A = complex<T>(n_acoefs[origin_coef], 0);
                B = complex<T>(n_acoefs[origin_coef+1], 0);
                C = complex<T>(n_acoefs[origin_coef+2], 0);

                auto first_comp_quad = quadr(A, B, C);
                D = complex<T>(first_comp_quad.first);
                E = complex<T>(first_comp_quad.second);

                /**< Make required substitutions, factorization again */
                complex<T> mul_tmp(un_bandwith, 0);
                complex<T> num_tmp(1, 0);
                complex<T> fr_tmp, result;
                A = complex<T>(1, 0);
                fr_tmp = num_tmp / D;
                fr_tmp = -fr_tmp;
                result = fr_tmp * mul_tmp;
                B = result;
                C = complex<T>(un_centrfreq * un_centrfreq, 0);

                auto second_comp_quad = quadr(A, B, C);
                D = complex<T>(second_comp_quad.first);
                E = complex<T>(second_comp_quad.second);

                /**< Determine final values for new coefficients */
                un_acoefs[new_coef] = 1.;
                un_acoefs[new_coef+1] = -2. * D.getReal();
                complex<T> fc_real = D * D.conj();
                un_acoefs[new_coef+2] = fc_real.getReal();
                un_acoefs[new_coef+3] = 1.;
                un_acoefs[new_coef+4] = -2. * E.getReal();
                complex<T> sc_real = E * E.conj();
                un_acoefs[new_coef+5] = sc_real.getReal();
            }

            /*<  */
            A = complex<T>(n_bcoefs[origin_coef], 0);
            B = complex<T>(n_bcoefs[origin_coef+1], 0);
            C = complex<T>(n_bcoefs[origin_coef+2], 0);

            auto first_comp_quad = quadr(A, B, C);
            D = std::move(first_comp_quad.first);
            E = std::move(first_comp_quad.second);

            /**< Make required substitutions, factorization again */
            complex<T> mul_tmp(un_bandwith, 0);
            complex<T> num_tmp(1, 0);
            complex<T> fr_tmp, result;
            A = complex<T>(1, 0);
            fr_tmp = num_tmp / D;
            fr_tmp = -fr_tmp;
            result = fr_tmp * mul_tmp;
            B = result;
            C = complex<T>(un_centrfreq * un_centrfreq, 0);

            un_bcoefs[new_coef] = 1.;
            un_bcoefs[new_coef+1] = -2. * D.getReal();
            complex<T> fc_real = std::move(D * D.conj());
            un_bcoefs[new_coef+2] = fc_real.getReal();
            un_bcoefs[new_coef+3] = 1.;
            un_bcoefs[new_coef+4] = -2. * E.getReal();
            complex<T> sc_real = std::move(E * E.conj());
            un_bcoefs[new_coef+5] = sc_real.getReal();
        }
    }
}

/**
 * @brief
 */
template<typename T>
void CalcFilterCoefs<T>::BPCoefsUnnorm(T un_bandwith, T un_centrfreq)
{
    std::size_t origin_qd_count,                         /**< Original number of quads values */
                origin_order;                            /**< Original order */
    std::size_t origin_coef, new_coef, pos_start;        /**< Counters */
    std::size_t size_coef;                               /**< Size vector value */
    complex<T> A, B, C, D, E;                            /**< Temp complex value */

    /**
      Store the original number of the order,
      new order will be twice large to original
    */
    origin_order = m_order;
    origin_qd_count = (origin_order + 1) / 2;
    m_order = origin_order * 2;
    /**<  */
    size_coef = 3*origin_order;

    if(FillZeroCoeffs(un_acoefs, un_bcoefs, size_coef))
    {
        /**< If original order is odd convert first order factor to quadratic,
        * pos_start indicate start point for loop */
        if(origin_order % 2)
        {
            un_acoefs[0] = n_acoefs[1];
            un_acoefs[1] = un_bandwith * n_acoefs[2];
            un_acoefs[2] = n_acoefs[1] * un_centrfreq * un_centrfreq;
            un_bcoefs[0] = n_bcoefs[1];
            un_bcoefs[1] = un_bandwith * n_bcoefs[2];
            un_bcoefs[2] = n_bcoefs[1] * un_centrfreq * un_centrfreq;
            pos_start = 1;
        }
        else
        {
            pos_start = 0;
        }

        /**< */
        for(std::size_t qd_count = pos_start; qd_count < origin_qd_count; qd_count++)
        {
            /**< origin_coef - for indexing original coefficients, every 3 per original coeffs from origin_qd_count
             *
             */

            origin_coef = qd_count * 3;
            new_coef = qd_count * 6 - pos_start * 3;

            if(n_acoefs[origin_coef] == 0)
            {
                un_acoefs[new_coef] = 0.;
                un_acoefs[new_coef+1] = std::sqrt(n_acoefs[origin_coef+2]) * un_bandwith;
                un_acoefs[new_coef+2] = 0.;
                un_acoefs[new_coef+3] = 0.;
                un_acoefs[new_coef+4] = std::sqrt(n_acoefs[origin_coef+2]) * un_bandwith;
                un_acoefs[new_coef+5] = 0.;
            }
            /**< */
            else
            {
                /**< Convert coefficients to complex, then factorization */
                A = complex<T>(n_acoefs[origin_coef], 0);
                B = complex<T>(n_acoefs[origin_coef+1], 0);
                C = complex<T>(n_acoefs[origin_coef+2], 0);

                auto first_comp_quad = quadr(A, B, C);
                D = complex<T>(first_comp_quad.first);
                E = complex<T>(first_comp_quad.second);

                /**< Make required substitutions, factorization again */
                complex<T> mul_tmp(un_bandwith, 0);
                complex<T> num_tmp(1, 0);
                complex<T> fr_tmp, result;
                A = complex<T>(1, 0);
                fr_tmp = num_tmp / D;
                fr_tmp = -fr_tmp;
                result = fr_tmp * mul_tmp;
                B = result;
                C = complex<T>(un_centrfreq * un_centrfreq, 0);

                auto second_comp_quad = quadr(A, B, C);
                D = complex<T>(second_comp_quad.first);
                E = complex<T>(second_comp_quad.second);

                /**< Determine final values for new coefficients */
                un_acoefs[new_coef] = 1.;
                un_acoefs[new_coef+1] = -2. * D.getReal();
                complex<T> fc_real = D * D.conj();
                un_acoefs[new_coef+2] = fc_real.getReal();
                un_acoefs[new_coef+3] = 1.;
                un_acoefs[new_coef+4] = -2. * E.getReal();
                complex<T> sc_real = E * E.conj();
                un_acoefs[new_coef+5] = sc_real.getReal();
            }

            /**<  */
            A = complex<T>(n_bcoefs[origin_coef], 0);
            B = complex<T>(n_bcoefs[origin_coef+1], 0);
            C = complex<T>(n_bcoefs[origin_coef+2], 0);

            auto first_comp_quad = quadr(A, B, C);
            D = std::move(first_comp_quad.first);
            E = std::move(first_comp_quad.second);

            /**< Make required substitutions, factorization again */
            complex<T> mul_tmp(un_bandwith, 0);
            complex<T> num_tmp(1, 0);
            complex<T> fr_tmp, result;
            A = complex<T>(1, 0);
            fr_tmp = num_tmp / D;
            fr_tmp = -fr_tmp;
            result = fr_tmp * mul_tmp;
            B = result;
            C = complex<T>(un_centrfreq * un_centrfreq, 0);

            un_bcoefs[new_coef] = 1.;
            un_bcoefs[new_coef+1] = -2. * D.getReal();
            complex<T> fc_real = std::move(D * D.conj());
            un_bcoefs[new_coef+2] = fc_real.getReal();
            un_bcoefs[new_coef+3] = 1.;
            un_bcoefs[new_coef+4] = -2. * E.getReal();
            complex<T> sc_real = std::move(E * E.conj());
            un_bcoefs[new_coef+5] = sc_real.getReal();
        }
    }
}

/**
 * @brief Method the denormalization of coefficients by frequency and
 *        impedance for high-pass filter type:
 *        \f[
 *             S(h) = \frac{\omega_0}{s} \bigg( \frac{\omega_{PB}}{s_n} \bigg)
 *           ]\f
 *        The gain constant multiplied by \f$(A_2/B_2)\f$
 */
template<typename T>
void CalcFilterCoefs<T>::HPCoefsUnnorm(T freq)
{
    /**< Number of coefficient, and position start */
    std::size_t coef_numb,
            pos_start,
            qd_count;

    if(FillZeroCoeffs(un_acoefs, un_bcoefs, m_order))
    {
        /**< First order type, if odd, set start position */
        if(m_order % 2)
        {
            m_fparam.gain *= (n_acoefs[2]/n_bcoefs[2]);
            un_acoefs[2] = freq*n_acoefs[1]/n_acoefs[2];
            un_acoefs[1] = 1.;
            un_bcoefs[2] = freq/n_bcoefs[2];
            pos_start = 1;
        }
        else
        {
            pos_start = 0;
        }

        for(qd_count = pos_start; qd_count < (m_order+1)/2; qd_count++)
        {
            coef_numb = qd_count*3;
            m_fparam.gain *= (n_acoefs[2]/n_bcoefs[2]);
            un_acoefs[coef_numb+1] = n_acoefs[coef_numb+1] * (freq/n_acoefs[coef_numb+2]);
            un_acoefs[coef_numb+2] = freq * freq * n_acoefs[coef_numb]/n_acoefs[coef_numb+2];
            un_acoefs[coef_numb] = 1.;
            un_bcoefs[coef_numb+1] = n_bcoefs[coef_numb+1] * (freq/n_bcoefs[coef_numb+2]);
            un_bcoefs[coef_numb+2] = freq * freq * n_bcoefs[coef_numb]/n_bcoefs[coef_numb+2];
            un_bcoefs[coef_numb] = 1.;
        }
    }
}

/**
 * @brief Method the denormalization of coefficients by frequency and
 *        impedance for low-pass filter type:
 *        \f[
 *             Z_s = z_0 \times Z_n(\omega_{PB} \times  s_n)
 *            ]\f
 *        and denormalized pole and zero defined as:
 *        \f[
 *             s = p \times z(\omega_{PB} \times p_n, \omega_{PB} \times z_n)
 *           ]\f
 *        The gain constant is unchanged.
 *        See ECE 6414: Continuous Time Filters(P. Allen)
 */
template<typename T>
void CalcFilterCoefs<T>::LPCoefsUnnorm(T freq)
{
    std::size_t nb_coef,
            qd_count,
            ps_start;

    if(FillZeroCoeffs(un_acoefs, un_bcoefs, m_order))
    {
        /**< First check order type, if odd, set start position to 1 else to 0 */
        if(m_order % 2)
        {
            un_acoefs[2] = n_acoefs[2]*freq;
            un_bcoefs[2] = n_bcoefs[2]*freq;
            ps_start = 1;
        }
        else
        {
            ps_start = 0;
        }

        for(qd_count=ps_start; qd_count < (m_order+1)/2; ps_start++)
        {
            nb_coef = qd_count*3;
            un_acoefs[nb_coef+1] = n_acoefs[nb_coef+1]*freq;
            un_acoefs[nb_coef+2] = n_acoefs[nb_coef+2]*(freq*freq);
            un_bcoefs[nb_coef+1] = n_bcoefs[nb_coef+1]*freq;
            un_bcoefs[nb_coef+2] = n_bcoefs[nb_coef+2]*freq;
        }
    }
}

template<typename T=double>
class CalcAnalogCoefs: public CalcFilterCoefs<T>
{};

template<typename T=double>
class CalcDigIIRCoefs: public CalcFilterCoefs<T>
{};

}
#endif //GENFILTER_H
