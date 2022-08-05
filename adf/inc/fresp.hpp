// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef FRESPONSE_H
#define FRESPONSE_H
#include <cmath>
#include <vector>

#include "base.hpp"

namespace adf {

/*!
 * \brief The Scale enum -
 */
enum class Scale
{
    LOG, /*!< Logarithmic scale */
    LIN  /*!< Linear scale */
};

/*!
 * \brief The FResolution enum -
 */
enum class FResolution
{
    HIGH=1,
    MEDIUM=2,
    LOW=3
};

/*!
 * \class Response -
 */
template<typename T=double>
class Response
{
    //!Common fiels initialize
    T m_start_freq = 0.;    /*!< start frequency for response */
    T m_stop_freq = 0.;     /*!< end frequency for response */
    T m_gain_min = 0.;      /*!< min gain to display */
    std::size_t m_decades=0;
    std::size_t m_dec_pts=0;
    std::size_t m_tot_pts=0;
    Scale m_frq;
    Scale m_mag;
    
    //! Field for estimation during initialization and fields used for subsequent calculations
    std::size_t m_resolution=0;
    T m_ratio=0.;
    std::vector<T> m_freq;       /*!< frequency values */
    std::vector<T> m_magn;       /*!< output magnitude values */
    std::vector<T> m_angl;       /*!< output phase values */
    std::vector<T> m_edge_magn;  /*!< edge frequency magnitudes */
    std::vector<T> m_edge_angl;  /*!< edge frequency phases */
    std::size_t m_num_edge_freq;

    void ResponseInit();

public:
    Response(const FResolution& resolution,
             T start_freq=ADF_FREQ_MIN,
             T stop_freq=ADF_FREQ_MAX,
             T gain_min=ADF_GAIN_STOP)
        :m_start_freq(start_freq)
        ,m_stop_freq(stop_freq)
        ,m_gain_min(gain_min)
    {
        
        switch(resolution)
        {
        case FResolution::HIGH:
            m_resolution=1;
            break;
        case FResolution::MEDIUM:
            m_resolution=2;
            break;
        case FResolution::LOW:
            m_resolution=3;
            break;
        }

        ResponseInit();
    }

    void VectorFill();

    void respAnalog(std::vector<T> &a_coeff,
                    std::vector<T> &b_coeff,
                    const std::size_t order,
                    const T gain);
    //void respDigitalIIR();
    //void respDigitalFIR();

    std::vector<T> getFrequency(){return m_freq;}
    std::vector<T> getMagnitude(){return m_magn;}
    std::vector<T> getPhase(){return m_angl;}
    //std::vector<T> getEdgeMagnitude();
    //std::vector<T> getEdgePhase();
    //std::size_t getNumEdgeFrequency();
};

/*!
 * \brief ResponseInit -
 */
template<typename T>
void Response<T>::ResponseInit()
{
    //! Determine the numbers point of frequencies, depending on the resolution of the response
    m_tot_pts = static_cast<std::size_t>(ADF_MAX_PTR / static_cast<std::size_t>(std::pow(2., m_resolution-1)));
    //! Response min gain
    m_gain_min = 10.0 * std::floor(m_gain_min/10.01);
    //! Determine if magnitude should be linear or logarithmic
    if(m_gain_min >= -0.001)
        m_mag = Scale::LIN;
    else
        m_mag = Scale::LOG;
    //! Determine lin or log frequency scale
    m_ratio = static_cast<std::size_t>(m_stop_freq/m_start_freq);

    if(m_ratio >= 10)
    {
        m_frq = Scale::LOG;

        while(m_ratio > 2)
        {
            ++m_decades;
            m_ratio /= 10;
        }

        m_dec_pts = m_tot_pts / m_decades;
        m_tot_pts = m_decades * m_dec_pts + 1;
        m_stop_freq = m_start_freq * std::pow(10, m_decades);
    }
    else{
        m_frq = Scale::LIN;
        m_dec_pts = 0;
    }
}

/*!
 * \brief VectorFill -
 */
template<typename T>
void Response<T>::VectorFill()
{

    T delta;
    m_freq.reserve(m_tot_pts);
    m_magn.reserve(m_tot_pts);
    m_angl.reserve(m_tot_pts);

    for(auto ind = 0; ind<m_tot_pts; ++ind)
    {
        m_magn.push_back(0.);
        m_angl.push_back(0.);
    }

    if(m_frq == Scale::LIN)
    {
        delta = (m_stop_freq - m_start_freq) / (m_tot_pts - 1);

        m_freq.pop_back(m_start_freq);
        for(auto ind = 1; ind<m_tot_pts; ++ind)
        {
            m_freq.push_back(m_freq[ind-1] + delta);
        }
    }
    else
    {
        delta = std::pow(10, 1.0/m_dec_pts);

        m_freq.pop_back(m_start_freq);
        for(auto ind = 1; ind<m_tot_pts; ++ind)
        {
            m_freq.push_back(m_freq[ind-1] * delta);
        }
    }
}

/*!
 * \brief respAnalog -
 * \param a_coeff -
 * \param b_coeff -
 * \param order -
 * \param gain -
 */
template<typename T>
void Response<T>::respAnalog(std::vector<T>& a_coeff, std::vector<T>& b_coeff, const std::size_t order, const T gain)
{
    VectorFill();

    for(auto find=0; find<m_tot_pts; ++find)
    {
        m_magn[find] = gain;

        auto omega = ADF_PI_2 * m_magn[find];
        auto pow2omega = omega * omega;

        for(std::size_t qind=0; qind<(order+1)/2; ++qind)
        {
            auto cindx = qind * 3;
            //! Numerator
            auto real = a_coeff[cindx + 2] - a_coeff[cindx] * pow2omega;
            auto imag = a_coeff[cindx + 1] * omega;
            auto mag = std::sqrt(real*real + imag*imag);
            m_magn[find] *= mag;

            if(mag > 0)
            {
                m_angl[find] += std::atan2(imag, real);
            }
            //! Denominator
            real = b_coeff[cindx + 2] - b_coeff[cindx] * pow2omega;
            imag = b_coeff[cindx + 1] * omega;
            mag = std::sqrt(real*real + imag*imag);
            m_magn[find] /= mag;

            if(mag > 0)
            {
                m_angl[find] -= std::atan2(imag, real);
            }
        }

        m_angl[find] *= ADF_RAD2DEG;
    }

    if(m_mag == Scale::LOG)
    {
        for(auto find=0; find<m_tot_pts; ++find)
        {
            if(m_magn[find] < ADF_ZERO)
            {
                m_magn[find] = ADF_ZERO;
            }
            m_magn[find] = 20 * std::log10(m_magn[find]);
        }
    }
}

}

#endif //FRESPONSE_H
