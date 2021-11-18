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

template<typename T=double>
class Response
{
    /*Common fiels initialize*/
    T m_start_freq = 0.; //start frequency for response
    T m_stop_freq = 0.; //end frequency for response
    T m_gain_min = 0.; //min gain to display
    std::size_t m_decades=0;
    std::size_t m_dec_pts=0;
    std::size_t m_tot_pts=0;
    Scale m_frq;
    Scale m_mag;
    
    /*Field for estimation during initialization and fields used for subsequent calculations*/
    std::size_t m_resolution=0;
    T m_ratio=0.;
    std::vector<T> freq; //frequency values
    std::vector<T> magn; //output magnitude values
    std::vector<T> angl; //output phase values
    std::vector<T> edge_magn; //edge frequency magnitudes
    std::vector<T> edge_angl; //edge frequency phases
    std::size_t num_edge_freq;

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
    }

    void VectorFill();

    void respAnalog();
    void respDigitalIIR();
    void respDigitalFIR();

    std::vector<T> getFrequency();
    std::vector<T> getMagnitude();
    std::vector<T> getPhase();
    //std::vector<T> getEdgeMagnitude();
    //std::vector<T> getEdgePhase();
    //std::size_t getNumEdgeFrequency();
};

template<typename T>
void Response<T>::ResponseInit()
{
    //Determine the numbers point of frequencies, depending on the resolution of the response
    m_tot_pts = static_cast<std::size_t>(ADF_MAX_PTR / static_cast<std::size_t>(std::pow(2., m_resolution-1)));
    //Response min gain
    m_gain_min = 10.0 * std::floor(m_gain_min/10.01);
    //Determine if magnitude should be lin or log
    if(m_gain_min >= -0.001)
        m_mag = Scale::LIN;
    else
        m_mag = Scale::LOG;
    //Determine lin or log frequency scale
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

template<typename T>
void Response<T>::VectorFill()
{

    T delta;
    freq.reserve(m_tot_pts);
    magn.reserve(m_tot_pts);
    angl.reserve(m_tot_pts);

    for(auto ind = 0; ind<m_tot_pts; ++ind)
    {
        magn.push_back(0.);
        angl.push_back(0.);
    }

    if(m_frq == Scale::LIN)
    {
        delta = (m_stop_freq - m_start_freq) / (m_tot_pts - 1);

        freq.pop_back(m_start_freq);
        for(auto ind = 1; ind<m_tot_pts; ++ind)
        {
            freq.push_back(freq[ind-1] + delta);
        }
    }
    else
    {
        delta = std::pow(10, 1.0/m_dec_pts);

        freq.pop_back(m_start_freq);
        for(auto ind = 1; ind<m_tot_pts; ++ind)
        {
            freq.push_back(freq[ind-1] * delta);
        }
    }
}

}

#endif //FRESPONSE_H
