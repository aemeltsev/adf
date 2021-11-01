#ifndef ADF_TESTS_H
#define ADF_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <gtest/gtest.h>

#include "inc/adfilter.hpp"

TEST(adfBttrwTst, BLPOrder)
{
    //Arrange
    constexpr int32_t blp_order = 3;
    adf::FiltParam<double> blp_fparam;

    blp_fparam.freq_passband.first = 0.1591549431;
    blp_fparam.freq_passband.second = 0.0;
    blp_fparam.freq_stopband.first = 0.3183098862;
    blp_fparam.freq_stopband.second = 0.0;
    blp_fparam.gain_passband.first = -1.0;
    blp_fparam.gain_passband.second = 0.0;
    blp_fparam.gain_stopband.first = -12.0;
    blp_fparam.gain_stopband.second = 0.0;

    adf::FilterType blp_filter = adf::FilterType::LPF;
    adf::ApproxType blp_apprx = adf::ApproxType::BUTTER;
    adf::CalcFilterCoefs<double> blp_coefs(blp_fparam, blp_filter, blp_apprx);

    //Act

    //Assert
    ASSERT_EQ(blp_coefs.getFilterOrder(), blp_order);
}

#endif //ADF_TESTS_H
