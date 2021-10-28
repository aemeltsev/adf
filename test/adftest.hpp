#ifndef ADF_TESTS_H
#define ADF_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <gtest/gtest.h>

#include "adfilter.hpp"

namespace adf_testing
{

class adfBttrwTst{

    adf::FiltParam m_fparam;
    adf::FilterSelect m_filter;
    adf::ApproxSelect m_apprx;
    adf::CalcFilterCoefs<> m_coefs;

public:
    adfBttrwTst();

};
}

#endif //ADF_TESTS_H
