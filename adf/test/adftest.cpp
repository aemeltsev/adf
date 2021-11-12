#ifndef ADF_TESTS_H
#define ADF_TESTS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cmath>
#include <gtest/gtest.h>

#include "inc/anfilter.hpp"

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
    auto m_coefs = new adf::AnalogFilter<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->setOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproximationType(), m_apprx);
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
    auto m_coefs = new adf::AnalogFilter<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->setOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproximationType(), m_apprx);
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
    auto m_coefs = new adf::AnalogFilter<double>(m_fparam, m_ftype, m_apprx);

    //Act
    m_coefs->setOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), m_ftype);
    ASSERT_EQ(m_coefs->getApproximationType(), m_apprx);
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
    auto m_coefs = new adf::AnalogFilter<double>(fparam, ftype, fapprx);

    //Act
    m_coefs->setOrder();

    //Assert
    ASSERT_EQ(m_coefs->getFilterType(), ftype);
    ASSERT_EQ(m_coefs->getApproximationType(), fapprx);
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

    auto afilter = new adf::AnalogFilter<>(fparam, ftype, fapprx);

    //Act
    afilter->setOrder();
    afilter->NormalizeCoefficients();
    auto&& ncoeffs = afilter->getNormalizeCoefficients();

    //Assert
    ASSERT_EQ(afilter->getFilterOrder(), order);
    ASSERT_TRUE(ncoeffs.first.size() != 0);
    ASSERT_TRUE(ncoeffs.second.size() != 0);
    ASSERT_EQ(ncoeffs.first.size(), cfsize);
    ASSERT_EQ(ncoeffs.second.size(), cfsize);
    delete afilter;
}

TEST(BttrwCoeffFill, BNcoeffCheck)
{
    constexpr std::size_t order = 5;
    constexpr int32_t cfsize = 9;
    std::vector<double> avec{0.0, 0.0, 1.23412,
                             0.0, 0.0, 1.52305,
                             0.0, 0.0, 1.52305};
    std::vector<double> bvec{0.0, 1.00000, 1.23412,
                             1.00000, 0.76273, 1.52305,
                             1.00000, 1.9968, 1.52305};

    auto ftype = adf::FilterType::LPF;
    auto fapprx = adf::ApproxType::BUTTER;
    adf::FiltParam<> fparam;
    tstFilterConf(ftype, fparam);
    auto afilter = new adf::AnalogFilter<>(fparam, ftype, fapprx);

    //Act
    afilter->setOrder();
    afilter->NormalizeCoefficients();
    auto&& ncoeffs = afilter->getNormalizeCoefficients();

    //Assert
    ASSERT_EQ(afilter->getFilterOrder(), order);
    ASSERT_TRUE(ncoeffs.first.size() > 0);
    ASSERT_TRUE(ncoeffs.second.size() > 0);
    ASSERT_EQ(ncoeffs.first.size(), cfsize);
    ASSERT_EQ(ncoeffs.second.size(), cfsize);
    for(auto i=0; i<cfsize; ++i)
    {
        ASSERT_EQ(std::round(ncoeffs.first[i]*10000)/10000, std::round(avec[i]*10000)/10000);
        ASSERT_EQ(std::round(ncoeffs.second[i]*10000)/10000, std::round(bvec[i]*10000)/10000);
    }
    delete afilter;
}

TEST(CHebyUnCoeff, CHebyUnCoeff)
{
    //Arrange
    constexpr int32_t order = 4;
    constexpr int32_t cfsize = 7;

    std::vector<double> navec{0.0, 0.0, 0.95046,
                              0.0, 0.0, 0.24336,
                              0.0};
    std::vector<double> nbvec{1.0, 0.23826, 0.95046,
                              1.0, 0.57521, 0.24335,
                              0.0};

    std::vector<double> unavec{0.0, 0.0, 1.23412,
                             0.0, 0.0, 1.52305,
                             0.0, 0.0, 1.52305};
    std::vector<double> unbvec{0.0, 1.00000, 1.23412,
                             1.00000, 0.76273, 1.52305,
                             1.00000, 1.9968, 1.52305};

    adf::FiltParam<double> fparam;
    auto ftype = adf::FilterType::HPF;
    auto fapprx = adf::ApproxType::CHEBY;

    fparam.freq_passband.first = 2000.0;//wp1
    fparam.freq_passband.second = 0.0;
    fparam.freq_stopband.first = 800.0;//ws1
    fparam.freq_stopband.second = 0.0;
    fparam.gain_passband.first = -1.5;
    fparam.gain_passband.second = 0.0;
    fparam.gain_stopband.first = -40.0;
    fparam.gain_stopband.second = 0.0;

    auto afilter = new adf::AnalogFilter<double>(fparam, ftype, fapprx);

    //Act
    afilter->setOrder();
    afilter->NormalizeCoefficients();
    auto&& ncoeffs = afilter->getNormalizeCoefficients();
    afilter->DenormalizeCoefficients();
    auto&& uncoeffs = afilter->getDenormalizeCoefficients();

    //Assert
    ASSERT_EQ(afilter->getFilterOrder(), order);
    ASSERT_TRUE(ncoeffs.first.size() > 0);
    ASSERT_TRUE(ncoeffs.second.size() > 0);
    ASSERT_EQ(ncoeffs.first.size(), cfsize);
    ASSERT_EQ(ncoeffs.second.size(), cfsize);
    for(auto i=0; i<cfsize; ++i)
    {
        ASSERT_EQ(std::round(ncoeffs.first[i]*10000)/10000, std::round(navec[i]*10000)/10000);
        ASSERT_EQ(std::round(ncoeffs.second[i]*10000)/10000, std::round(nbvec[i]*10000)/10000);
    }
    ASSERT_EQ(uncoeffs.first.size(), cfsize);
    ASSERT_EQ(uncoeffs.second.size(), cfsize);
    delete afilter;
}

TEST(EllUnCoeff, EllUnCoeffCheck)
{
    //Arrange
    constexpr int32_t order = 4;
    constexpr int32_t cfsize = 7;

    adf::FiltParam<double> fparam;
    auto ftype = adf::FilterType::LPF;
    auto fapprx = adf::ApproxType::BUTTER;

    fparam.freq_passband.first = 8620.0;//wp1
    fparam.freq_passband.second = 0.0;
    fparam.freq_stopband.first = 16880.0;//ws1
    fparam.freq_stopband.second = 0.0;
    fparam.gain_passband.first = -1.0;
    fparam.gain_passband.second = 0.0;
    fparam.gain_stopband.first = -17.0;
    fparam.gain_stopband.second = 0.0;

    auto afilter = new adf::AnalogFilter<double>(fparam, ftype, fapprx);

    //Act
    afilter->setOrder();
    afilter->NormalizeCoefficients();
    auto&& ncoeffs = afilter->getNormalizeCoefficients();
    afilter->DenormalizeCoefficients();
    auto&& uncoeffs = afilter->getDenormalizeCoefficients();

    //Assert
    ASSERT_EQ(afilter->getFilterOrder(), order);
}

#endif //ADF_TESTS_H
