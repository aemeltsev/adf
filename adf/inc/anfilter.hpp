// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef ANALOGFILTER_H
#define ANALOGFILTER_H

#include "genfilter.hpp"

namespace adf {
/*!
 * \class AnalogFilter - analog active filter solver class
 */
template<typename T=double>
class AnalogFilter
{
    CalcCoeffs<T>* m_calccoeffs;
    FiltParam<T> m_fparam;
    FilterType m_ftype;
    ApproxType m_atype;
    std::size_t m_order = 0;
    std::vector<T> n_acoefs, n_bcoefs; /*!< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /*!< to unnormalise coefs */

    /*!
     * \brief setFilterParam Filling data fields of the
     * \param g_passband - the passband gain ripple
     * \param g_stopband - the stopband gain ripple
     * \param f_passband - the passband edge frequency
     * \param f_stopband - the stopband edge frequency
     * \param fsamp - sample frequency
     * \param gain - gain miltiplier
     */
    void setFilterParam(FiltParam<T>& other)
    {

        m_fparam.gain_passband = std::move(other.gain_passband);
        m_fparam.gain_stopband = std::move(other.gain_stopband);
        m_fparam.freq_passband = std::move(other.freq_passband);
        m_fparam.freq_stopband = std::move(other.freq_stopband);
        m_fparam.fsamp = other.fsamp;
        m_fparam.gain = other.gain;
    }

    /*!
     * \brief FillZeroCoeffs - Filling normalized vectors with default values
     * \param avec input reference to vector of "a" coefficients for filling
     * \param bvec input reference to vector of "b" coefficients for filling
     * \param num the order
     * \return bool variable, success if order value not out of range
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
    /*!
     * \brief AnalogFilter<T> - ctor
     * \param fparam - input specification for filter parameters
     * \param ftype - type filter select
     * \param atype - select the approximation methods
     */
    explicit AnalogFilter<T>(FiltParam<T>& fparam, FilterType& ftype, ApproxType& atype)
        :m_ftype(ftype)
        ,m_atype(atype)
    {
        setFilterParam(fparam);
        m_calccoeffs = new CalcCoeffs<T>(fparam, ftype, atype);
    }

    /*!
      * \brief dtor
      */
    ~AnalogFilter<T>()
    {
        delete m_calccoeffs;
    }

    /*!
     * \brief setOrder - add filter order
     */
    void setOrder()
    {
        m_calccoeffs->FilterOrder();
        m_order = m_calccoeffs->getFilterOrder();
    }

    /*!
     * \brief NormalizeCoefficients - Calculation and filling the vectors of coefficients depending on the approximation method
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
            throw std::range_error(ADF_ERROR("The order value to out of range"));
        }
    }

    /*!
     * \brief  - Denormalization consists in the transition
     *           from normalized parameters to true ones.
     */
    void DenormalizeCoefficients()
    {
        T freq, /*!< Stored value for unnormaliztion frequency */
          BW,   /*!< For transformation to bandstop or bandpass type */
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

    /*!
     * \brief getFilterOrder
     * \return filter order value
     */
    std::size_t getFilterOrder() const
    {
        return m_order;
    }

    /*!
     * \brief getApproximationType
     * \return return approximation type value(integer)
     */
    ApproxType getApproximationType()
    {
        return m_atype;
    }

    /*!
     * \brief getFilterType
     * \return return filter type value(integer)
     */
    FilterType getFilterType()
    {
        return m_ftype;
    }

    /*!
     * \brief getNormalizeCoefficients
     * \return pair of vectors with a, b normalized coefficients
     */
    std::pair<std::vector<T>, std::vector<T>> getNormalizeCoefficients()
    {
        return std::make_pair(n_acoefs, n_bcoefs);
    }

    /*!
     * \brief getDenormalizeCoefficients
     * \return pair of vectors with a, b unnormalized coefficients
     */
    std::pair<std::vector<T>, std::vector<T>> getDenormalizeCoefficients()
    {
        return std::make_pair(un_acoefs, un_bcoefs);
    }
};

} //adf

#endif //ANALOGFILTER_H
