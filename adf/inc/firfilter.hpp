// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef DIGFIRFILTER_H
#define DIGFIRFILTER_H

#include "genfilter.hpp"

namespace adf {

enum class WindowType
{
    RectangularWindow,
    HammingWindow,
    HanningWindow,
    BlackmanWindow,
    KaiserWindow,
    ParksMcClellanWindow,
    BartlettWindow,
    Undefined
};

template<typename T>
class DigFIRFilter
{
private:

    FiltParam<T> m_fparam;
    FilterType m_ftype;
    ApproxType m_atype;
    std::size_t m_order = 0;
    std::vector<T> n_acoefs, n_bcoefs; /*!< to normalise coefs */
    //std::vector<T> un_acoefs, un_bcoefs; /*!< to unnormalise coefs */

    T CalcFilterLen(const WindowType &w_type);
    void CalcIdealFIRCoeffs();
public:

    DigFIRFilter(FilterType& f_type, FiltParam<T>& f_param, ApproxType& atype)
    {

    }

    void CalcDigFIRCoeffs(const WindowType& w_type, std::size_t order = 0);

    void RectangularWindowCoeffs();
    void HammingWindowCoeffs();
    void HanningWindowCoeffs();
    void BlackmanWindowCoeffs();
    void KaiserWindowCoeffs();
    void BartlettWindowCoeffs();
    void ParksMcClellanWindowCoeffs();
    void setParksMcClellanParameters();
    void RemezInterchange();
    void ComputeLagrange();
    void EstimateFreqResponse();
    void ComputeCoeffs();
    void MultWinIdealCoeffs();
};

/*!
 * \brief
 * \param
 * \result
 */
template<typename T>
T DigFIRFilter<T>::CalcFilterLen(const WindowType& w_type)
{
    T result = 0.;
    T freq_delta, freq_delta_lower, freq_delta_upper;
    if(m_ftype == FilterType::LPF || m_ftype == FilterType::HPF)
    {
        freq_delta = std::fabs(m_fparam.freq_stopband.first()
                               - m_fparam.freq_passband.first()) / m_fparam.fsamp;
    }
    else if(m_ftype == FilterType::PBF || m_ftype == FilterType::PBF)
    {
        freq_delta_lower = std::fabs(m_fparam.freq_stopband.first()
                                     - m_fparam.freq_passband.first()) / m_fparam.fsamp;
        freq_delta_upper = std::fabs(m_fparam.freq_stopband.second()
                                     - m_fparam.freq_passband.second()) / m_fparam.fsamp;
        if(freq_delta_lower > freq_delta_upper)
        {
            freq_delta = freq_delta_upper;
        }
        else{
            freq_delta = freq_delta_lower;
        }
    }

    auto sb_error = std::pow(10.0, 0,05 * m_fparam.gain_stopband.first());
    auto pb_error = 1 - std::pow(10.0, 0,05 * m_fparam.gain_passband.first());

    T min_error;
    if(w_type == WindowType::RectangularWindow ||
       w_type == WindowType::HammingWindow ||
       w_type == WindowType::HanningWindow ||
       w_type == WindowType::BlackmanWindow ||
       w_type == WindowType::BartlettWindow ||
       w_type == WindowType::KaiserWindow)
    {
        if(sb_error < pb_error){
            min_error = sb_error;
        } else{
            min_error = pb_error;
        }

        auto db_error = -20 * std::log10(min_error);
        if(db_error > 50){
            result = 0.1102 * (db_error - 8.7);
        }
        else if(db_error > 21){
            result = 0.5842 * std::pow((db_error - 21), 0.4) + 0.07886 * (db_error - 21);
        }
        else{
            result = 0.;
        }

        if(db_error > 21){
            m_order = static_cast<std::size_t>(std::ceil((db_error - 7.95) / (2.285 * freq_delta))
                                               + 1);
        }
        else{
            m_order = static_cast<std::size_t>(std::ceil((5.794 / freq_delta) + 1));
        }

        if(!(m_order % 2)){
            ++m_order;
        }

    } else if(w_type == WindowType::ParksMcClellanWindow){
        sb_error = std::log10(sb_error);
        pb_error = std::log10(pb_error);
        freq_delta /= ADF_PI_2;

        auto k1 = (0.005309 * pb_error * pb_error + 0.07114 * pb_error - 0.4761) * sb_error
                - (0.00266 * pb_error * pb_error + 0.5941 * pb_error + 0.4278);
        auto k2 = 0.51244 * (pb_error - sb_error) = 11.012;
        m_order = static_cast<std::size_t>(std::ceil(((k1 - k2 * freq_delta * freq_delta)
                                                      / freq_delta) + 1.5));

        if(!(m_order % 2)){
            ++m_order;
        }
        result = 0.;
    }
    return result;
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::CalcIdealFIRCoeffs()
{
    auto tau = static_cast<T>(m_order  - 1) / 2;

    auto cut_freq1 = (m_fparam.freq_stopband.first() + m_fparam.freq_passband.first()) / (2 * m_fparam.fsamp);
    auto cut_freq2 = (m_fparam.freq_stopband.second() + m_fparam.freq_passband.second()) / (2 * m_fparam.fsamp);

    switch(m_ftype)
    {
    case FilterType::LPF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = cut_freq1 / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = std::sin(cut_freq1 * t) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::HPF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (ADF_PI - cut_freq1) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(ADF_PI * t) - std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::PBF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (cut_freq2 - cut_freq1) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(cut_freq2 * t) - std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::SBF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (ADF_PI + cut_freq1 - cut_freq2) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(ADF_PI * t) - std::sin(cut_freq2 * t)
                               + std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;
    }
}

template<typename T>
void DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType& w_type, std::size_t order)
{
    auto beta = CalcFilterLen(w_type);

    if(order < 0) return;
    if(order > 0)
    {
        m_order = order;
    }

    n_acoefs.reserve(m_order);
    n_bcoefs.reserve(m_order);
    n_acoefs.assign(m_order, 0);
    n_bcoefs.assign(m_order, 0);

    m_fparam.gain = 1.0;

    if(w_type != WindowType::ParksMcClellanWindow)
    {
        CalcIdealFIRCoeffs();
    }

    switch(w_type)
    {
    case WindowType::RectangularWindow :
        break;

    case WindowType::BartlettWindow :
        break;

    case WindowType::BlackmanWindow :
        break;

    case WindowType::HammingWindow :
        break;

    case WindowType::HanningWindow :
        break;

    case WindowType::KaiserWindow :
        break;

    case WindowType::ParksMcClellanWindow :
        break;

    case WindowType::Undefined :
        ADF_ERROR("Error in DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType&, std::size_t): The unknown window type");
        return;

    default:
        ADF_ERROR("Error in DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType&, std::size_t): The window type value, is empty");
        return;
    }

    if(w_type != WindowType::ParksMcClellanWindow)
    {
        MultWinIdealCoeffs();
    }
}

} //namespace adf

#endif //DIGFIRFILTER_H