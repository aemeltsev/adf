/**
  * License text...
*/
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
    explicit complex(T re, T im=0) noexcept
        :m_re(re)
        ,m_im(im)
    {}
    //copy ctor
    complex(const complex& other)
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
};

/**
 *  @brief - square root of complex number
 *  @param - unique pointer to complex value
 *  @return - returns the square root of complex number
 */
template<typename T>
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
template<typename T>
inline std::pair<complex<T>, complex<T>> quadr(
        complex<T>& p_a,
        complex<T>& p_b,
        complex<T>& p_c)
{
    complex<T> a2_var,
               ac4_tmp, ac4_var,
               sq1_tmp, sq2_tmp, sq_var,
               fr_root, sc_root,
               pb_tmp = p_b;
    complex<T> fr_const(2., 0.);
    complex<T> sc_const(4., 0.);

    a2_var = fr_const * p_a; /**< 2*a */

    ac4_tmp = sc_const * p_a;
    ac4_var = ac4_tmp * p_c; /**< 4*a*c */

    sq1_tmp = p_b * p_b;
    sq2_tmp = sq1_tmp - ac4_var;
    sq_var = adf::sqrt(sq2_tmp); /**< sqrt(b*b - 4*a*c */

    pb_tmp = -pb_tmp;
    fr_root = pb_tmp + sq_var; /**< for first root */
    sc_root = pb_tmp - sq_var; /**< for second root */
    return std::make_pair(fr_root/a2_var, sc_root/a2_var);
}

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
}
#endif // COMPLEX_H
