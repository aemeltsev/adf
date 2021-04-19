/**
  * License text...
*/
#ifndef COMPLEX_H
#define COMPLEX_H
#include <iostream>
#include <cstdint>
#include <memory>
#include <cmath>
#include <utility>

namespace adf {

template<class T>
class complex
{
    /**
     * @brief Working variable for saving re, im
     */
    T m_re, m_im;

public:
    //ctors
    complex() noexcept {m_re=m_im=0;}
    complex(T re, T im=0)
        :m_re(re)
        ,m_im(im)
    {}
    //copy ctor
    explicit complex(const complex& other)
        :m_re(other.m_re)
        ,m_im(other.m_im)
    {}
    //copy assign
    complex<T> &operator=(const complex<T>& other)
    {
        if(this != other){
            m_re = other.m_re;
            m_im = other.m_im;
        }
        return *this;
    }
    //move ctor
    complex(const complex&& other) noexcept
        :m_re(std::move(other.m_re))
        ,m_im(std::move(other.m_im))
    {}
    //move assign
    complex<T> &operator=(const complex<T>&& other)
    {
        if(this != other){
            m_re = std::move(other.m_re);
            m_im = std::move(other.m_im);
        }
        return *this;
    }

    /**
     * @brief getters
     */
    T getReal() const {return m_re;}
    T getImag() const {return m_im;}

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
    template<typename U> friend bool operator==(complex<U> &p_fr, complex<U> &p_sc);
    template<typename U> friend bool operator!=(complex<U> &p_fr, complex<U> &p_sc);

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

    complex<T> operator+() const {return *this;}
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
    friend std::ofstream &operator<<(std::ostream &p_out, const complex<T> &p_val)
    {
        return p_out << "(" << p_val.getReal() << "," << p_val.getImag() << ")";
    }
};

/**
 *  @brief - square root of complex number
 *  @param - unique pointer to complex value
 *  @return - returns the square root of complex number
 */
template<class T>
inline complex<T> sqrt(complex<T>& p_val)
{
    T real = std::sqrt(p_val.mag()) * std::cos(p_val.arg()/2.);
    T imag = std::sqrt(p_val.mag()) * std::sin(p_val.arg()/2.);
    return complex<T>(real, imag);
}

/**
 *  @brief - factors quadratic equation with cmplx coefficients
 *  @param - a*x^2
 *  @param - b*x
 *  @param - c
 *  @return - std::pairs with complex roots
 */
template<class T>
inline std::pair<complex<T>, complex<T>> quadr(
        complex<T>& p_a,
        complex<T>& p_b,
        complex<T>& p_c)
{
    complex<T> a2_var, ac4_var, sq_var;
    complex<T> fr_const(2., 0.);
    complex<T> sc_const(4., 0.);
    a2_var = fr_const * p_a;
    ac4_var = sc_const * p_a * p_c;
    sq_var = adf::sqrt(p_b * p_b - ac4_var);
    return std::make_pair(((-p_b) + sq_var)/a2_var, ((-p_b) - sq_var)/a2_var);
}

#include "complex.cpp"
}
#endif // COMPLEX_H
