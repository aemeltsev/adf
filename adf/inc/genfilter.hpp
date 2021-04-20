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
 * @brief The FilterSelect enum - the enumeration class for the type filter select
 */
enum class FilterSelect
{
    LPF,
    HPF,
    PBF,
    SBF
};

/**
 * @brief The ApproxSelect enum - the enumeration class for select the approximation methods
 */
enum class ApproxSelect
{
    BUTTER,
    CHEBY,
    ICHEBY,
    ELLIPT
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
    std::unique_ptr<FiltParam<T>> m_fparam;
    FilterSelect m_sfilter;
    ApproxSelect m_sapprox;
    uint16_t m_order;
    T m_gain;
    T CommonKernel();
    T FreqNorm();
    void FilterOrder();

protected:
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */
    void setFiltParam(
            std::pair<T, T>& g_passband, /**< The pasband gain ripple */
            std::pair<T, T>& g_stopband, /**< The stopband gain ripple */
            std::pair<T, T>& f_passband,
            std::pair<T, T>& f_stopband,
            T fsamp,
            T gain
            /*,int16_t order*/);
    void setTypeFilter(FilterSelect& sfilter);
    void setApproxFilter(ApproxSelect& sapprox);

    void ButterApprox();
    void ChebyApprox();
    void ElliptApprox();
    void IChebyApprox();
    void NormalCoefs();
    void BSCoefsUnnorm(T un_bandwith, T un_centrfreq);
    void BPCoefsUnnorm(T un_bandwith, T un_centrfreq);
    void HPCoefsUnnorm(T freq);
    void LPCoefsUnnorm(T freq);
    void UnnormCoefs();

public:
    explicit CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam, FilterSelect &fselect, ApproxSelect &sapprox) noexcept;
    CalcFilterCoefs() noexcept;
};

#include "genfilter.cpp"
}
#endif //GENFILTER_H
