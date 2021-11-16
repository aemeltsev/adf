#ifndef COMMON_H
#define COMMON_H

#include <cstdint>
#include <string>

namespace adf {

/* constants */
constexpr double ADF_PI(3.14159265358979);   /**< pi */
constexpr double ADF_PI_2(6.28318530717959); /**< 2 * pi */
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

constexpr long MAX_TERMS(100);
constexpr double ERR_SMALL(1e-15);

constexpr double ADF_RAD2DEG(180.0/ADF_PI);

constexpr long double ADF_ZERO(1e-30);
constexpr long ADF_MAX_PTR(640);

/* errors */
std::string error(const std::string& err, const char* func, const char* file, int16_t line)
{
    std::string result = err;
    result.append(1, ' ').append(func);
    result.append(1, ' ').append(file);
    result.append(1, ' ').append(std::to_string(line));
    return result;
}

#define ADF_ERROR(msg) error(msg, __FUNCTION__, __FILE__, __LINE__)

}

#endif //COMMON_H
