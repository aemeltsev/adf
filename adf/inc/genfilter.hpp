#ifndef GENFILTER_H
#define GENFILTER_H
#include <cstdint>
#include <memory>
#include <cmath>
#include <vector>

#include "advmath.hpp"
#include "complex.hpp"
#include "errornum.hpp"

namespace adf {

enum class FilterSelect
{
    LPF,
    HPF,
    PBF,
    SBF
};

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
  int16_t order;
};

/**
 * @brief
 */
template<class T>
class CalcFilterCoefs
{
private:
    std::unique_ptr<FiltParam<T>> m_fparam;
    FilterSelect m_sfilter;
    ApproxSelect m_sapprox;

protected:
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */
    void setFiltParam(
            std::pair<T, T>& g_passband,
            std::pair<T, T>& g_stopband,
            std::pair<T, T>& f_passband,
            std::pair<T, T>& f_stopband,
            T fsamp,
            T gain,
            int16_t order);
    void setTypeFilter(FilterSelect& sfilter);
    void setApproxFilter(ApproxSelect& sapprox);
    T CommonKernel();

public:
    explicit CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam, FilterSelect &fselect, ApproxSelect &sapprox) noexcept;
    CalcFilterCoefs() noexcept;

    void FilterOrder();
    void NormalCoefs();
    void ButterApprox();
    void ChebyApprox();
    void ElliptApprox();
    void IChebyApprox();
    void UnnormCoefs();
    void BSCoefsUnnorm();
    void BPCoefsUnnorm();
    void HPCoefsUnnorm();
    void LPCoefsUnnorm();
};

#include "genfilter.cpp"
}
#endif //GENFILTER_H
