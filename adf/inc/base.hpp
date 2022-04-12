// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef COMMON_H
#define COMMON_H

#include <iostream>
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
/* errors */

template <class T> T zero(T t)
{
    T result;
    if(std::is_same<decltype(t), char>::value)
    {
        result = static_cast<char>(0);
    }
    else if(std::is_same<decltype(t), short>::value)
    {
        result = static_cast<short>(0);
    }
    else if(std::is_same<decltype(t), int>::value)
    {
        result = static_cast<int>(0);
    }
    else if(std::is_same<decltype(t), long>::value)
    {
        result = static_cast<long>(0);
    }
    else if(std::is_same<decltype(t), float>::value)
    {
        result = static_cast<float>(0.0);
    }
    else if(std::is_same<decltype(t), double>::value)
    {
        result = static_cast<double>(0.0);
    }
    else{
        result = T() - T();
    }
    return result;
}

template<class T> T one(T t)
{
    T result;
    if(std::is_same<decltype(t), char>::value)
    {
        result = static_cast<char>(1);
    }
    else if(std::is_same<decltype(t), short>::value)
    {
        result = static_cast<short>(1);
    }
    else if(std::is_same<decltype(t), int>::value)
    {
        result = static_cast<int>(1);
    }
    else if(std::is_same<decltype(t), long>::value)
    {
        result = static_cast<long>(1);
    }
    else if(std::is_same<decltype(t), float>::value)
    {
        result = static_cast<float>(1.0);
    }
    else if(std::is_same<decltype(t), double>::value)
    {
        result = static_cast<double>(1.0);
    }
    else{
        std::cerr << " one() not implemented " << std::endl;
    }
    return result;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

}

#endif //COMMON_H
