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

/* constants */
constexpr double ADF_PI(3.1415926535898);   /**< pi */
constexpr double ADF_PI_2(6.2831853071796); /**< 2 * pi */
constexpr double ADF_IOTA(1.0E-3);          /**< small increment */
constexpr double ADF_FREQ_MIN(1.0E-3);      /**< min freq */
constexpr double ADF_FREQ_MAX(1.0E12);      /**< max freq */
constexpr double ADF_GAIN_PASS(-1.0E-2);    /**< min passband gain */
constexpr double ADF_GAIN_TRAN(-3.0103);    /**< min sb, max pb gain */
constexpr double ADF_GAIN_STOP(-2.0E02);    /**< max stopband gain */
constexpr double ADF_GRID_MULT(16);         /**< starting grid mult */
constexpr double ADF_MAX_BANDS(10);         /**< max bands for PM */
constexpr double ADF_SYM_ODD(1);            /**< type 1 FIR */
constexpr double ADF_SYM_EVEN(2);           /**< type 2 FIR */
constexpr double ADF_ASYM_ODD(3);           /**< type 3 FIR */
constexpr double ADF_ASYM_EVEN(4);          /**< type 4 FIR */

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

protected:
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */

public:
    explicit CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam) noexcept;
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
