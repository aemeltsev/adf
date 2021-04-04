#include "complex.hpp"
namespace adf {

template<class T>
inline bool operator==(complex<T> &p_fr, complex<T> &p_sc)
{
    return ((p_fr.getReal()==p_sc.getReal())&&(p_fr.getImag()==p_sc.getImag())) ? true : false;
}

template<class T>
inline bool operator!=(complex<T> &p_fr, complex<T> &p_sc)
{
    return (p_fr==p_sc) ? true : false;
}

template<class T>
inline complex<T> operator+(complex<T> &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr.m_re+p_sc.m_re), (p_fr.m_im+p_sc.m_im));
}

template<class T>
inline complex<T> operator+(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr+p_sc.m_re), p_sc.m_im);
}

template<class T>
inline complex<T> operator+(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re+p_fr), p_sc.m_im);
}

template<class T>
inline complex<T> operator-(complex<T> &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr.m_re+p_sc.m_re), (p_fr.m_im+p_sc.m_im));
}

template<class T>
inline complex<T> operator-(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr-p_sc.m_re), p_sc.m_im);
}

template<class T>
inline complex<T> operator-(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re-p_fr), p_sc.m_im);
}

template<class T>
inline complex<T> operator*(complex<T> &p_fr, complex<T> &p_sc)
{
    T real = p_fr.getReal()*p_sc.getReal() - p_fr.getImag()*p_sc.getImag();
    T imag = p_fr.getReal()*p_sc.getImag() + p_sc.getReal()*p_fr.getImag();
    return complex<T>(real, imag);
}

template<class T>
inline complex<T> operator*(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_sc.m_re*p_fr), (p_sc.m_im*p_fr));
}

template<class T>
inline complex<T> operator*(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re*p_fr), (p_sc.m_im*p_fr));
}

template<class T>
inline complex<T> operator/(complex<T> &p_fr, complex<T> &p_sc)
{
    T real = (p_fr.getReal()*p_sc.getImag() + p_fr.getImag()*p_sc.getImag())/(p_sc.getReal()*p_sc.getReal() + p_sc.getImag()*p_sc.getImag());
    T imag = (p_sc.getReal()*p_fr.getImag() - p_fr.getReal()*p_sc.getImag())/(p_sc.getReal()*p_sc.getReal() + p_sc.getImag()*p_sc.getImag());
    return complex<T>(real, imag);
}

template<class T>
complex<T> &complex<T>::operator+=(complex<T> &p_val)
{
    m_re += p_val.m_re;
    m_im += p_val.m_im;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator+=(T p_val)
{
    m_re += p_val;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator-=(complex<T> &p_val)
{
    m_re -= p_val.m_re;
    m_im -= p_val.m_im;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator-=(T p_val)
{
    m_re -= p_val;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator*=(complex<T> &p_val)
{
    m_re = m_re*p_val.getReal() - m_im*p_val.getImag();
    m_im = m_re*p_val.getImag() + p_val.getReal()*m_im;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator*=(T p_val)
{
    m_re *= p_val;
    m_im *= p_val;
    return *this;
}

template<class T>
complex<T> &complex<T>::operator/=(complex<T> &p_val)
{
    m_re = (m_re*p_val.getReal() + m_im*p_val.getImag())/(p_val.getReal()*p_val.getReal() + p_val.getImag()*p_val.getImag());
    m_im = (p_val.getReal()*m_im - m_re*p_val.getImag())/(p_val.getReal()*p_val.getReal() + p_val.getImag()*p_val.getImag());
    return *this;
}

template<class T>
complex<T> &complex<T>::operator/=(T p_val)
{
    m_re /= p_val;
    m_im /= p_val;
    return *this;
}

}
