// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef COMPLEX_H
#define COMPLEX_H
#include <iostream>
#include <cmath>

namespace adf {

template<typename T>
class complex
{
    /**
     * @brief Working variable for store re, im
     */
    T m_re, m_im;

public:
    //ctors
    complex() noexcept {m_re=m_im=0;}
    complex(T re, T im=0) noexcept
        :m_re(re)
        ,m_im(im)
    {}
    //copy ctor
    complex(const complex<T>& other)
        :m_re(other.m_re)
        ,m_im(other.m_im)
    {}
    //copy assign
    complex<T> &operator=(const complex<T>& other)
    {
        if(this != &other){
            m_re = other.m_re;
            m_im = other.m_im;
        }
        return (*this);
    }
    //move ctor
    complex(const complex&& other) noexcept
        :m_re(std::move(other.m_re))
        ,m_im(std::move(other.m_im))
    {}
    //move assign
    complex<T> &operator=(complex<T>&& other)
    {
        if(this != &other){
            m_re = std::move(other.m_re);
            m_im = std::move(other.m_im);
        }
        return (*this);
    }

    /**
     * @brief getters
     */
    T getReal() const {return m_re;}
    T getImag() const {return m_im;}
    void setReal(T val) {m_re = val;}
    void setImag(T val) {m_im = val;}

    //math methods
    //conj() - returns complex conj of complex number
    complex<T> conj(){return complex<T>(m_re,-m_im);}
    T norm() const {return (m_re*m_re+m_im*m_im);}
    //arg() - returns angle (radians) of complex number
    T arg() const {return (!m_re&&!m_im) ? 0 : std::atan2(m_im, m_re);}
    //mag() - returns the magnitude of complex number
    T mag() const {return std::sqrt(m_re*m_re + m_im*m_im);}

    /**
     * @brief overloadings
     */
    template<typename U> friend complex<U> operator+(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator+(U &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator+(complex<U> &p_sc, U &p_fr);
    template<typename U> friend complex<U> operator-(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator-(U &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator-(complex<U> &p_sc, U &p_fr);
    template<typename U> friend complex<U> operator*(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator*(U &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator*(complex<U> &p_sc, U &p_fr);
    template<typename U> friend complex<U> operator/(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator/(U &p_fr, complex<U> &p_sc);
    template<typename U> friend complex<U> operator/(complex<U> &p_sc, U &p_fr);

    template<typename U> friend bool operator==(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend bool operator!=(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend bool operator==(const complex<U> &p_fr, const complex<U> &p_sc);
    template<typename U> friend bool operator!=(const complex<U> &p_fr, const complex<U> &p_sc);

    complex<T> operator+() {return *this;}
    complex<T> operator+() const {return *this;}
    complex<T> operator-() {return complex<T>(-m_re, -m_im);}
    complex<T> operator-() const {return complex<T>(-m_re, -m_im);}

    complex<T> &operator+=(complex<T> &p_val);
    complex<T> &operator+=(T p_val);
    complex<T> &operator-=(complex<T> &p_val);
    complex<T> &operator-=(T p_val);
    complex<T> &operator*=(complex<T> &p_val);
    complex<T> &operator*=(T p_val);
    complex<T> &operator/=(complex<T> &p_val);
    complex<T> &operator/=(T p_val);
    /**
     * @brief out stream
     */
    template<typename U>
    friend std::ostream &operator<<(std::ostream &p_out, const complex<U> &p_val);
    template<T> friend class FFT;
};

template<typename T>
std::ostream &operator<<(std::ostream &p_out, const complex<T> &p_val)
{
    p_out << p_val.getReal();
    T im = p_val.getImag();
    if(im < 0)
        p_out << im << "i";
    else if(im > 0)
        p_out << "+" << im << "i";
    return p_out;
}

template<typename T>
inline bool operator==(complex<T> &p_fr, complex<T> &p_sc)
{
    return ((p_fr.getReal()==p_sc.getReal())&&(p_fr.getImag()==p_sc.getImag())) ? true : false;
}

template<typename T>
inline bool operator!=(complex<T> &p_fr, complex<T> &p_sc)
{
    return (p_fr==p_sc) ? true : false;
}

template<typename T>
inline bool operator==(const complex<T> &p_fr, const complex<T> &p_sc)
{
    return ((p_fr.getReal()==p_sc.getReal())&&(p_fr.getImag()==p_sc.getImag())) ? true : false;
}

template<typename T>
inline bool operator!=(const complex<T> &p_fr, const complex<T> &p_sc)
{
    return (p_fr==p_sc) ? true : false;
}

template<typename T>
inline complex<T> operator+(complex<T> &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr.m_re+p_sc.m_re), (p_fr.m_im+p_sc.m_im));
}

template<typename T>
inline complex<T> operator+(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr+p_sc.m_re), p_sc.m_im);
}

template<typename T>
inline complex<T> operator+(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re+p_fr), p_sc.m_im);
}

template<typename T>
inline complex<T> operator-(complex<T> &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr.m_re+p_sc.m_re), (p_fr.m_im+p_sc.m_im));
}

template<typename T>
inline complex<T> operator-(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_fr-p_sc.m_re), p_sc.m_im);
}

template<typename T>
inline complex<T> operator-(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re-p_fr), p_sc.m_im);
}

template<typename T>
inline complex<T> operator*(complex<T> &p_fr, complex<T> &p_sc)
{
    T real = p_fr.getReal()*p_sc.getReal() - p_fr.getImag()*p_sc.getImag();
    T imag = p_fr.getReal()*p_sc.getImag() + p_sc.getReal()*p_fr.getImag();
    return complex<T>(real, imag);
}

template<typename T>
inline complex<T> operator*(T &p_fr, complex<T> &p_sc)
{
    return complex<T>((p_sc.m_re*p_fr), (p_sc.m_im*p_fr));
}

template<typename T>
inline complex<T> operator*(complex<T> &p_sc, T &p_fr)
{
    return complex<T>((p_sc.m_re*p_fr), (p_sc.m_im*p_fr));
}

template<typename T>
inline complex<T> operator/(complex<T> &p_fr, complex<T> &p_sc)
{
    T real = (p_fr.getReal()*p_sc.getImag() + p_fr.getImag()*p_sc.getImag())/(p_sc.getReal()*p_sc.getReal() + p_sc.getImag()*p_sc.getImag());
    T imag = (p_sc.getReal()*p_fr.getImag() - p_fr.getReal()*p_sc.getImag())/(p_sc.getReal()*p_sc.getReal() + p_sc.getImag()*p_sc.getImag());
    return complex<T>(real, imag);
}

template<typename T>
complex<T> &complex<T>::operator+=(complex<T> &p_val)
{
    m_re += p_val.m_re;
    m_im += p_val.m_im;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator+=(T p_val)
{
    m_re += p_val;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator-=(complex<T> &p_val)
{
    m_re -= p_val.m_re;
    m_im -= p_val.m_im;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator-=(T p_val)
{
    m_re -= p_val;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator*=(complex<T> &p_val)
{
    m_re = m_re*p_val.getReal() - m_im*p_val.getImag();
    m_im = m_re*p_val.getImag() + p_val.getReal()*m_im;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator*=(T p_val)
{
    m_re *= p_val;
    m_im *= p_val;
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator/=(complex<T> &p_val)
{
    m_re = (m_re*p_val.getReal() + m_im*p_val.getImag())/(p_val.getReal()*p_val.getReal() + p_val.getImag()*p_val.getImag());
    m_im = (p_val.getReal()*m_im - m_re*p_val.getImag())/(p_val.getReal()*p_val.getReal() + p_val.getImag()*p_val.getImag());
    return *this;
}

template<typename T>
complex<T> &complex<T>::operator/=(T p_val)
{
    m_re /= p_val;
    m_im /= p_val;
    return *this;
}

} //namespace adf
#endif // COMPLEX_H
