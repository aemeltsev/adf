#ifndef ADF_TESTS_H
#define ADF_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>
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
    constexpr int32_t order = 5;

    adf::FiltParam<double> m_fparam;
    auto m_ftype = adf::FilterType::LPF;
    auto m_apprx = adf::ApproxType::BUTTER;
    tstFilterConf(m_ftype, m_fparam);
    auto m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), order);
    delete m_coefs;
}

TEST(CHebyTstSuite, CHHPApprx)
{
    //Arrange
    constexpr int32_t order = 4;

    adf::FiltParam<double> m_fparam;
    auto m_ftype = adf::FilterType::HPF;
    auto m_apprx = adf::ApproxType::CHEBY;
    tstFilterConf(m_ftype, m_fparam);
    auto m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), order);
    delete m_coefs;
}

TEST(ICHebyTstSuite, ICHBPApprx)
{
    //Arrange
    constexpr int32_t order = 3;

    adf::FiltParam<double> m_fparam;
    auto m_ftype = adf::FilterType::PBF;
    auto m_apprx = adf::ApproxType::ICHEBY;
    tstFilterConf(m_ftype, m_fparam);
    auto m_coefs = new adf::CalcFilterCoefs<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproxType(), m_apprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), order);
    delete m_coefs;
}

TEST(ElliptTstSuite, ELLBSApprx)
{
    //Arrange
    constexpr int32_t order = 3;

    adf::FiltParam<double> fparam;
    auto ftype = adf::FilterType::SBF;
    auto fapprx = adf::ApproxType::ELLIPT;
    tstFilterConf(ftype, fparam);
    auto m_coefs = new adf::CalcFilterCoefs<double>(fparam, ftype, fapprx);

    //Act
    m_coefs->makeFilterOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), ftype);
    ASSERT_EQ(m_coefs->getApproxType(), fapprx);
    ASSERT_EQ(m_coefs->getFilterOrder(), order);
    delete m_coefs;
}

TEST(BttrwCoeffFill, BNcoeffSize)
{
    constexpr std::size_t order = 7;
    constexpr std::size_t cfsize = 12;

    auto ftype = adf::FilterType::LPF;
    auto fapprx = adf::ApproxType::BUTTER;
    adf::FiltParam<> fparam;

    fparam.freq_passband.first = 12000.0;//wp1
    fparam.freq_passband.second = 0.0;
    fparam.freq_stopband.first = 20000.0;//ws1
    fparam.freq_stopband.second = 0.0;
    fparam.gain_passband.first = -0.5;
    fparam.gain_passband.second = 0.0;
    fparam.gain_stopband.first = -21.0;
    fparam.gain_stopband.second = 0.0;

    auto fcoefs = new adf::CalcFilterCoefs<>(fparam, ftype, fapprx);

    //Act
    fcoefs->makeFilterOrder();
    fcoefs->makeNormalCoefs();

    //Assert
    ASSERT_EQ(fcoefs->getFilterOrder(), order);
    ASSERT_TRUE(fcoefs->normACoefsSize() != 0);
    ASSERT_EQ(fcoefs->normACoefsSize(), cfsize);
    delete fcoefs;
}

TEST(BttrwCoeffFill, BNcoeffCheck)
{
    constexpr std::size_t order = 5;
    constexpr int32_t cfsize = 9;
    std::vector<double> avec{0.0, 0.0, 1.23412,
                             0.0, 0.0, 1.52305,
                             0.0, 0.0, 1.52305};
    std::vector<double> bvec{0.0, 1.000, 1.234,
                             1.000, 1.997, 1.523,
                             1.000, 0.7627, 1.523};

    auto ftype = adf::FilterType::LPF;
    auto fapprx = adf::ApproxType::BUTTER;
    adf::FiltParam<> fparam;
    tstFilterConf(ftype, fparam);
    auto fcoefs = new adf::CalcFilterCoefs<>(fparam, ftype, fapprx);

    //Act
    fcoefs->makeFilterOrder();
    fcoefs->makeNormalCoefs();
    auto acoef = fcoefs->normACoefs();

    //Assert
    for(auto i=0; i<cfsize; ++i)
    {
        ASSERT_EQ(std::round(acoef[i]*10000)/10000, std::round(avec[i]*10000)/10000);
    }
    ASSERT_EQ(fcoefs->getFilterOrder(), order);
    ASSERT_TRUE(fcoefs->normACoefsSize() > 0);
    ASSERT_EQ(fcoefs->normACoefsSize(), cfsize);
    delete fcoefs;
}
#endif //ADF_TESTS_H
