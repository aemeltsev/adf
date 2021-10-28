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
 * @brief Container for input filter params, used double type by default
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
    std::size_t m_order;
    std::size_t m_gain;
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */

    T CommonKernel();
    T FreqNorm();
    void FilterOrder();

protected:
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

    std::size_t getFilterOrder();
    ApproxType getApproxType();
    FilterType getFilterType();
    std::vector<T> normACoefs();
    std::vector<T> normBCoefs();
    std::vector<T> unnormACoefs();
    std::vector<T> unnormBCoefs();

public:
    explicit CalcFilterCoefs(const FiltParam<T> &fparam, const FilterType &fselect, const ApproxType &sapprox) noexcept;
    CalcFilterCoefs() noexcept;

    void setFiltParam(
            std::pair<T, T> &g_passband, /**< The pasband gain ripple */
            std::pair<T, T> &g_stopband, /**< The stopband gain ripple */
            std::pair<T, T> &f_passband,
            std::pair<T, T> &f_stopband,
            T fsamp,
            T gain
            /*,int16_t order*/);
    void setTypeFilter(FilterType& sfilter);
    void setApproxFilter(ApproxType& sapprox);
};

template<typename T=double>
class CalcAnalogCoefs: public CalcFilterCoefs<T>
{

};

template<typename T=double>
class CalcDigIIRCoefs: public CalcFilterCoefs<T>
{

};


}
#endif //GENFILTER_H
