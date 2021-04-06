#include "genfilter.hpp"
namespace adf {

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs(std::unique_ptr<FiltParam<T>> fparam) noexcept
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
