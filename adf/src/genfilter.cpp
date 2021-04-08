#include "genfilter.hpp"
namespace adf {

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam, FilterSelect &sfilter, ApproxSelect &sapprox) noexcept
    :m_sfilter(sfilter)
    ,m_sapprox(sapprox)
{
    m_fparam = std::make_unique<FiltParam<T>>();
    m_fparam->gain_passband = std::move(fparam->gain_passband);
    m_fparam->gain_stopband = std::move(fparam->gain_stopband);
    m_fparam->freq_passband = std::move(fparam->freq_passband);
    m_fparam->freq_stopband = std::move(fparam->freq_stopband);
    m_fparam->fsamp = std::move(fparam->fsamp);
    m_fparam->gain = std::move(fparam->gain);
    m_fparam->order = std::move(fparam->order);

}

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs() noexcept
{
    m_fparam = std::make_unique<FiltParam<T>>();
}

template<typename T>
void CalcFilterCoefs<T>::setFiltParam(
        std::pair<T, T>& g_passband,
        std::pair<T, T>& g_stopband,
        std::pair<T, T>& f_passband,
        std::pair<T, T>& f_stopband,
        T fsamp,
        T gain,
        int16_t order)
{
    m_fparam->gain_passband = std::move(g_passband);
    m_fparam->gain_stopband = std::move(g_stopband);
    m_fparam->freq_passband = std::move(f_passband);
    m_fparam->freq_stopband = std::move(f_stopband);
    m_fparam->fsamp = std::move(fsamp);
    m_fparam->gain = std::move(gain);
    m_fparam->order = std::move(order);
}

template<typename T>
void CalcFilterCoefs<T>::setTypeFilter(FilterSelect& sfilter)
{
    m_sfilter = sfilter;
}

template<typename T>
void CalcFilterCoefs<T>::setApproxFilter(ApproxSelect& sapprox)
{
    m_sapprox = sapprox;
}

/**
 * @brief This the common value for all filter approximation methods
 *        \f$(\varepsilon_s / \varepsilon_p)
 *        \f$(\varepsilon_s = 10.0^{-0.1*a_s}-1) stopband gain adjustment factor and
 *        passband gain ratio \f$(\varepsilon_p = 10.0^{-0.1*a_p}-1)
 * @return Ratio value of the suppression \f$(R_s)dB/\f$(R_p)dB
 */
template<typename T>
T CalcFilterCoefs<T>::CommonKernel()
{
    return ((std::pow(10.0,-0.1*m_fparam->gain_stopband.first)-1)/
            (std::pow(10.0,-0.1*m_fparam->gain_passband.first)-1));
}

template<typename T>
void CalcFilterCoefs<T>::FilterOrder()
{
    T kernel, ratio,
      order,
      wp1, wp2, ws1, ws2; /**< edge frequency variables */

    wp1 = m_fparam->f_passband.first;
    wp2 = m_fparam->f_passband.second;
    ws1 = m_fparam->f_stopband.first;
    ws2 = m_fparam->f_stopband.second;

    switch(m_sfilter)
    {
    case FilterSelect::LPF:
        ratio = ws1/wp1;
        break;
    case FilterSelect::HPF:
        ratio = wp1/ws1;
        break;
    case FilterSelect::PBF:
        if(ws1 > (wp1 * wp2) / ws2)
        {
            ws2 = (wp1 * wp2) / ws1;
            m_fparam->f_stopband.first = ws2;
        }
        else
        {
            ws1 = (wp1 * wp2) / ws2;
            m_fparam->f_stopband.second = ws1;
        }
        ratio = (ws2 - ws1) / (wp2 - wp1);
        break;
    case FilterSelect::SBF:
        if(wp1 > (ws1 * ws2) / wp2)
        {
            wp2 = (ws1 * ws2) / wp1;
            m_fparam->f_passband.first = wp2;
        }
        else
        {
            wp1 = (ws1 * ws2) / wp2;
            m_fparam->f_passband.second = wp1;
        }
        ratio = (wp2 - wp1) / (ws2 - ws1);
        break;
    default: return /* error code */;
    }

    kernel = CommonKernel();

    switch(m_sapprox)
    {
    case ApproxSelect::BUTTER:
        order = std::log10(kernel)/(2 * std::log10(ratio));
        break;
    case ApproxSelect::CHEBY:
    case ApproxSelect::ICHEBY:
        order = acosh(std::sqrt(kernel))/acosh(ratio);
        break;
    case ApproxSelect::ELLIPT:
        T ratio_const = 1/ratio;
        if(ratio_const > .9999) return /* error code */;
        T kernel_const = 1/std::sqrt(kernel);
        if(ratio_const < 2e-8) return /* error code */;
    }
}

template<typename T>
void CalcFilterCoefs<T>::NormalCoefs()
{

}

template<typename T>
void CalcFilterCoefs<T>::ButterApprox()
{

}

template<typename T>
void CalcFilterCoefs<T>::ChebyApprox()
{

}

template<typename T>
void CalcFilterCoefs<T>::ElliptApprox()
{

}

template<typename T>
void CalcFilterCoefs<T>::IChebyApprox()
{

}

template<typename T>
void CalcFilterCoefs<T>::UnnormCoefs()
{

}

template<typename T>
void CalcFilterCoefs<T>::BSCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::BPCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::HPCoefsUnnorm()
{

}

template<typename T>
void CalcFilterCoefs<T>::LPCoefsUnnorm()
{

}

}
