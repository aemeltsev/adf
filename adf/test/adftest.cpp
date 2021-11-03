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
        m_fparam.freq_passband.first = 1000.0;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 2000.0;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -0.5;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -21.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::HPF:
        m_fparam.freq_passband.first = 2000.0;
        m_fparam.freq_passband.second = 0.0;
        m_fparam.freq_stopband.first = 800.0;
        m_fparam.freq_stopband.second = 0.0;
        m_fparam.gain_passband.first = -1.5;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -40.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::PBF:
        m_fparam.freq_passband.first = 100.0; //wp1
        m_fparam.freq_passband.second = 200.0; //wp2
        m_fparam.freq_stopband.first = 50.0; //ws1
        m_fparam.freq_stopband.second = 400.0; //ws2
        m_fparam.gain_passband.first = -0.5;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -33.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    case adf::FilterType::SBF:
        m_fparam.freq_passband.first = 50.0; //wp1
        m_fparam.freq_passband.second = 72.0; //wp2
        m_fparam.freq_stopband.first = 58.0; //ws1
        m_fparam.freq_stopband.second = 62.0; //ws2
        m_fparam.gain_passband.first = -0.3;
        m_fparam.gain_passband.second = 0.0;
        m_fparam.gain_stopband.first = -50.0;
        m_fparam.gain_stopband.second = 0.0;
        break;
    }
}


TEST(BttrwTstSuite, BLPApprx)
{
    //Arrange
    constexpr int32_t blp_order = 5;

    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::LPF;
    adf::ApproxType m_apprx = adf::ApproxType::BUTTER;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), blp_order);
    delete m_coefs;
}

TEST(CHebyTstSuite, CHHPApprx)
{
    //Arrange
    constexpr int32_t blp_order = 4;

    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::HPF;
    adf::ApproxType m_apprx = adf::ApproxType::CHEBY;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), blp_order);
    delete m_coefs;
}

TEST(ICHebyTstSuite, ICHBPApprx)
{
    //Arrange
    constexpr int32_t blp_order = 3;

    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::PBF;
    adf::ApproxType m_apprx = adf::ApproxType::ICHEBY;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), blp_order);
    delete m_coefs;
}

TEST(ElliptTstSuite, ELLBSApprx)
{
    //Arrange
    constexpr int32_t blp_order = 3;

    adf::FiltParam<double> m_fparam;
    adf::FilterType m_ftype = adf::FilterType::SBF;
    adf::ApproxType m_apprx = adf::ApproxType::ELLIPT;
    adf::CalcFilterCoefs<double>* m_coefs;
    tstFilterConf(m_ftype, m_fparam);
    m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), blp_order);
    delete m_coefs;
}

#endif //ADF_TESTS_H
