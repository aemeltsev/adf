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
    MEDIUM,
    LOW
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

    void ResponseInit();

public:
    CalcResponse(const FiltResponse<T>& fresp)
    {

    }

};


}

#endif //FRESPONSE_H
