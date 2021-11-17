#ifndef FRESPONSE_H
#define FRESPONSE_H
#include <cmath>
#include <vector>

#include "base.hpp"
namespace adf {

enum class Scale
{
    LOG, //Logarithmic scale
    LIN  //Linear scale
};

enum class FResolution
{
    HIGH=1,
    MEDIUM=2,
    LOW=3
};

/**
 * @brief
 * @param freq
 * @param magn
 * @param angl
 * @param edge_magn
 * @param edge_angl
 * @param start_freq
 * @param stop_freq
 * @param gain_min
 * @param gain_min
 * @param num_edge_freq;
 * @param decades;
 * @param dec_pts;
 * @param tot_pts;
 * @param frq;
 * @param mag;
 */
template<typename T=double>
struct FiltResponse
{
    std::vector<T> freq;
    std::vector<T> magn;
    std::vector<T> angl;
    std::vector<T> edge_magn;
    std::vector<T> edge_angl;
    std::size_t num_edge_freq;
};

template<typename T=double>
class CalcResponse
{
    FiltResponse<T> m_fresp;
    T m_start_freq = 0.;
    T m_stop_freq = 0.;
    T m_gain_min = 0.;
    std::size_t m_decades=0;
    std::size_t m_dec_pts=0;
    std::size_t m_tot_pts=0;
    Scale m_frq;
    Scale m_mag;
    
    std::size_t m_resolution=0;
    T m_ratio=0.;

    void ResponseInit();

public:
    CalcResponse(FiltResponse<T>& fresp,
                 const FResolution& resolution,
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
    }

};

template<typename T>
void CalcResponse<T>::ResponseInit()
{
    m_tot_pts = static_cast<std::size_t>(ADF_MAX_PTR / static_cast<std::size_t>(std::pow(2., m_resolution-1)));
    m_gain_min = 10.0 * std::floor(m_gain_min/10.01);

    if(m_gain_min >= -0.001)
        m_mag = Scale::LIN;
    else
        m_mag = Scale::LOG;
    
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

}

#endif //FRESPONSE_H
