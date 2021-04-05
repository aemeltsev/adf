#include "genfilter.hpp"
namespace adf {

template<typename T>
CalcFilterCoefs<T>::CalcFilterCoefs(FiltParam<T>& fparam) noexcept
{
    m_fparam = std::move(fparam);
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
