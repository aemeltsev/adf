// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef DIGIIRFILTER_H
#define DIGIIRFILTER_H

#include "genfilter.hpp"

namespace adf {

template<typename T=double>
class DigIIRFilter
{
    CalcCoeffs<T>* m_calccoeffs;
    FiltParam<T> m_fparam;
    FilterType m_ftype;
    ApproxType m_atype;
    std::size_t m_order = 0;
    std::vector<T> n_acoefs, n_bcoefs; /**< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /**< to unnormalise coefs */

    /**
     * @brief setFilterParam Filling data fields of the
     * @param g_passband - the passband gain ripple
     * @param g_stopband - the stopband gain ripple
     * @param f_passband - the passband edge frequency
     * @param f_stopband - the stopband edge frequency
     * @param fsamp - sample frequency
     * @param gain - gain miltiplier
     */
    void setFilterParam(FiltParam<T>& other)
    {

        m_fparam.gain_passband = other.gain_passband;
        m_fparam.gain_stopband = other.gain_stopband;
        m_fparam.freq_passband = other.freq_passband;
        m_fparam.freq_stopband = other.freq_stopband;
        m_fparam.fsamp = other.fsamp;
        m_fparam.gain = other.gain;
    }

    /**
     * @brief FillZeroCoeffs - Filling normalized vectors with default values
     * @param avec input reference to vector of "a" coefficients for filling
     * @param bvec input reference to vector of "b" coefficients for filling
     * @param num the order
     * @return bool variable, success if order value not out of range
     */
    bool FillZeroCoeffs(std::vector<T>& avec, std::vector<T>& bvec, std::size_t num)
    {
        bool success = false;
        if(num > 0 && num < MAX_TERMS)
        {
            success = true;
            std::size_t number_coeffs = (3 * (num + 1)) / 2;
            avec.reserve(number_coeffs);
            bvec.reserve(number_coeffs);

            for(std::size_t ind=0; ind<number_coeffs; ++ind)
            {
                avec.push_back(0);
                bvec.push_back(0);
            }
        }
        return success;
    }

public:
    explicit DigIIRFilter<T>(FiltParam<T>& fparam, FilterType& ftype, ApproxType& atype)
    {
        if(ftype == FilterType::UNDEF || ftype >= FilterType::ERR)
        {
            throw std::invalid_argument(ADF_ERROR("DigIIRFilter : Undefined or error filter type"));
        }

        if(atype == ApproxType::UNDEF || atype == ApproxType::ERR)
        {
            throw std::invalid_argument(ADF_ERROR("DigIIRFilter : Undefined or error type of approximation"));
        }
        m_ftype = ftype;
        m_atype = atype;
        setFilterParam(fparam);
        m_calccoeffs = new CalcCoeffs<T>(fparam, ftype, atype);
    }

    ~DigIIRFilter<T>()
    {
        delete m_calccoeffs;
    }

    /**
     * @brief WarpFreq
     */
    void WarpFreq()
    {
        if(m_ftype == FilterType::LPF || m_ftype == FilterType::HPF){

            m_fparam.freq_passband.first = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_passband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.first = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_stopband.first / (2 * m_fparam.fsamp));
        }
        else if(m_ftype == FilterType::PBF || m_ftype == FilterType::SBF){

            m_fparam.freq_passband.first = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_passband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_passband.second = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_passband.second / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.first = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_stopband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.second = 2 * m_fparam.fsamp *
                                           std::tan(m_fparam.freq_stopband.second / (2 * m_fparam.fsamp));
        }
    }

    /**
     * @brief setOrder
     */
    void setOrder()
    {
        m_calccoeffs->FilterOrder();
        m_order = m_calccoeffs->getFilterOrder();
    }

    /**
     * @brief NormalizeCoefficients - Calculation and filling the vectors of coefficients depending on the approximation method
     */
    void NormalizeCoefficients()
    {
        if(FillZeroCoeffs(n_acoefs, n_bcoefs, m_order))
        {
            switch (m_atype)
            {
            case ApproxType::BUTTER:
                m_calccoeffs->ButterApprox(n_acoefs, n_bcoefs);
                break;
            case ApproxType::CHEBY:
                m_calccoeffs->ChebyApprox(n_acoefs, n_bcoefs);
                break;
            case ApproxType::ELLIPT:
                m_calccoeffs->ElliptApprox(n_acoefs, n_bcoefs);
                break;
            case ApproxType::ICHEBY:
                m_calccoeffs->IChebyApprox(n_acoefs, n_bcoefs);
                break;
            }
        }
        else
        {
            throw std::range_error(ADF_ERROR("DigIIRFilter : The order value to out of range"));
        }
    }

    /**
     * @brief DenormalizeCoefficients
     */
    void DenormalizeCoefficients()
    {
        T freq, // Stored value for unnormaliztion frequency
          BW,   /* For transformation to bandstop or bandpass type */
          Wo;

        if(m_atype == ApproxType::BUTTER
           || m_atype == ApproxType::CHEBY
           || m_atype == ApproxType::ELLIPT)
        {
            freq = m_fparam.freq_passband.first;
            Wo = std::sqrt(m_fparam.freq_passband.first * m_fparam.freq_passband.second);
            BW = m_fparam.freq_passband.second - m_fparam.freq_passband.first;
        }
        else if(m_atype == ApproxType::ICHEBY)
        {
            freq = m_fparam.freq_stopband.first;
            Wo = std::sqrt(m_fparam.freq_stopband.first * m_fparam.freq_stopband.second);
            BW = m_fparam.freq_passband.second - m_fparam.freq_passband.first;
        }

        switch (m_ftype)
        {
        case FilterType::LPF:
            if(FillZeroCoeffs(un_acoefs, un_bcoefs, m_order))
            {
                m_calccoeffs->LPCoefsUnnorm(n_acoefs, n_bcoefs, un_acoefs, un_bcoefs, freq);
            }
            break;
        case FilterType::HPF:
            if(FillZeroCoeffs(un_acoefs, un_bcoefs, m_order))
            {
                m_calccoeffs->HPCoefsUnnorm(n_acoefs, n_bcoefs, un_acoefs, un_bcoefs, freq);
            }
            break;
        case FilterType::PBF:
            m_calccoeffs->BPCoefsUnnorm(n_acoefs, n_bcoefs, un_acoefs, un_bcoefs, BW, Wo);
            break;
        case FilterType::SBF:
            m_calccoeffs->BSCoefsUnnorm(n_acoefs, n_bcoefs, un_acoefs, un_bcoefs, BW, Wo);
            break;
        }
    }

    /**
     * @brief BilTransform
     */
    void BilTransform()
    {
        auto c = 2 * m_fparam.fsamp;
        auto c_sq = c * c;

        auto half_order = (m_order + 1) / 2;
        auto start = 0;
        T az0, az1, az2,
          bz0, bz1, bz2;

        if(m_order % 2)
        {
            az0 = n_acoefs[2] + n_acoefs[1] * c;
            az1 = n_acoefs[2] - n_acoefs[1] * c;
            bz0 = n_bcoefs[2] + n_bcoefs[1] * c;
            bz1 = n_bcoefs[2] - n_bcoefs[1] * c;
            n_acoefs[0] = 1.;
            n_acoefs[1] = az1 / az0;
            n_acoefs[2] = 0.;
            n_bcoefs[0] = 1.;
            n_bcoefs[1] = bz1 / bz0;
            n_bcoefs[2] = 0.;
            m_fparam.gain *= (az0 / bz0);
            start = 1;
        }

        for(auto i = start; i < half_order; ++i)
        {
            auto j = 3 * i;
            az0 = n_acoefs[j] * c_sq + n_acoefs[j + 1] * c + n_acoefs[j + 2];
            az1 = 2 * (n_acoefs[j + 2] - n_acoefs[j] * c_sq);
            az2 = n_acoefs[j] * c_sq - n_acoefs[j + 1] * c + n_acoefs[j + 2];
            bz0 = n_bcoefs[j] * c_sq + n_bcoefs[j + 1] * c + n_bcoefs[j + 2];
            bz1 = 2 * (n_bcoefs[j + 2] - n_bcoefs[j] * c_sq);
            bz2 = n_bcoefs[j] * c_sq - n_bcoefs[j + 1] * c + n_bcoefs[j + 2];
            n_acoefs[j] = 1.;
            n_acoefs[j + 1] = az1 / az0;
            n_acoefs[j + 2] = az2 / az0;
            n_bcoefs[j] = 1.;
            n_bcoefs[j + 1] = bz1 / bz0;
            n_bcoefs[j + 2] = bz2 / bz0;
            m_fparam.gain *= (az0 / bz0);
        }
    }

    /**
     * @brief unWarpFreq
     */
    void unWarpFreq()
    {
        if(m_ftype == FilterType::LPF || m_ftype == FilterType::HPF){

            m_fparam.freq_passband.first = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_passband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.first = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_stopband.first / (2 * m_fparam.fsamp));
        }
        else if(m_ftype == FilterType::PBF || m_ftype == FilterType::SBF){

            m_fparam.freq_passband.first = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_passband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_passband.second = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_passband.second / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.first = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_stopband.first / (2 * m_fparam.fsamp));
            m_fparam.freq_stopband.second = 2 * m_fparam.fsamp *
                                           std::atan(m_fparam.freq_stopband.second / (2 * m_fparam.fsamp));
        }
    }

};

} //adf

#endif //DIGIIRFILTER_H
