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


template<typename T>
void CalcFilterCoefs<T>::FilterOrder()
{

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
