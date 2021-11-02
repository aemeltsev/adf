#ifndef ADF_TESTS_H
#define ADF_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <gtest/gtest.h>

#include "inc/adfilter.hpp"

void tstFilterConf(const adf::FilterType& m_ftype, adf::FiltParam<double>& m_fparam)
{
    switch (m_ftype)
    {
    case adf::FilterType::LPF:
        m_fparam.freq_passband.first = 0.1591549431;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 0.3183098862;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -1.0;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -12.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::HPF:
        m_fparam.freq_passband.first = 0.1591549431;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 0.3183098862;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -1.0;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -12.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::PBF:
        m_fparam.freq_passband.first = 0.1591549431;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 0.3183098862;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -1.0;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -12.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::SBF:
        m_fparam.freq_passband.first = 0.1591549431;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 0.3183098862;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -1.0;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -12.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    }
}


TEST(BttrwTstCase, BLPApprx)
{
    //Arrange
    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::LPF;
    adf::ApproxType m_apprx = adf::ApproxType::BUTTER;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act

    //Assert
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
}

TEST(BttrwTstSuite, BLPType)
{
    //Arrange
    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::LPF;
    adf::ApproxType m_apprx = adf::ApproxType::BUTTER;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
}

TEST(BttrwTstSuite, BLPOrder)
{
    //Arrange
    constexpr int32_t blp_order = 3;

    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::LPF;
    adf::ApproxType m_apprx = adf::ApproxType::BUTTER;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterOrder(), blp_order);
}





#endif //ADF_TESTS_H
