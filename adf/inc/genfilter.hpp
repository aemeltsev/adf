#ifndef GENFILTER_H
#define GENFILTER_H
#include <cstdint>
#include <cmath>
#include <vector>

#include "advmath.hpp"
#include "complex.hpp"
#include "errornum.hpp"

namespace adf {

/* #DEFINES/CONSTANTS: */
#define ADF_PI          3.1415926535898 /* pi */
#define ADF_PI_2        6.2831853071796 /* 2 * pi */
#define ADF_IOTA        1.0E-3          /* small increment */
#define ADF_FREQ_MIN    1.0E-3          /* min freq */
#define ADF_FREQ_MAX    1.0E12          /* max freq */
#define ADF_GAIN_PASS   -1.0E-2         /* min passband gain */
#define ADF_GAIN_TRAN   -3.0103         /* min sb, max pb gain */
#define ADF_GAIN_STOP   -2.0E02         /* max stopband gain */
#define ADF_GRID_MULT   16              /* starting grid mult */
#define ADF_MAX_BANDS   10              /* max bands for PM */
#define ADF_SYM_ODD     1               /* type 1 FIR */
#define ADF_SYM_EVEN    2               /* type 2 FIR */
#define ADF_ASYM_ODD    3               /* type 3 FIR */
#define ADF_ASYM_EVEN   4               /* type 4 FIR */

template<class T>
class CalcFilterCoefs
{

};

template<class T>
class generalized_filter
{
public:
    generalized_filter()
    {}

private:
    std::pair<T, T> gain_passband; /* apass_one + apass_two */
    std::pair<T, T> gain_stopband; /* astop_one + astop_two */
    std::pair<T, T> freq_passband; /* wpass_one + wpass_two */
    std::pair<T, T> freq_stopband; /* wstop_one + wstop_two */
    std::vector<T> acoefs, bcoefs; /* ref's to coef*/

    T fsamp; /* samp freq for dig filt */
    T gain; /* gain multiplier */
    int8_t order; /* order, length of filter */
    char select, approx, implem; /* selectivity, approximation and implementation */

    void BilinearTransform();

    void CalcAnalogCoefs();
    void CalcBartWinCoefs();
    void CalcBlckWinCoefs();
    void CalcButterCoefs();
    void CalcChebyCoefs();
    void CalcDigFIRCoefs();
    void CalcDigIIRCoefs();
    void CalcElliptCoefs();
    void CalcFilterCoefs();
    void CalcFilterOrder();
    void CalcHammWinCoefs();
    void CalcHannWinCoefs();
    void CalcIChebyCoefs();
    void CalcIdealFIRCoefs();
    void CalcKaisWinCoefs();
    void CalcNormalCoefs();
    void CalcParkMcclCoefs();
    void CalcRectWinCoefs();
    void ComputeCoefs();


};

}

#endif //GENFILTER_H
