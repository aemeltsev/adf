// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef POLYNOM_H
#define POLYNOM_H

#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
#include <initializer_list>
#include "base.hpp"
#include "complex.hpp"

#define MR 8             /*!< number of different fractional values to break (rare) limit cycles */
#define MT 10            /*!< number of steps to break (rare) limit cycles */
#define MAX_IT (MT * MR) /*!< maximal number of iterations */

namespace adf
{

/*!
 * \class Polynomial - parameterized class for polynomial working
 *        not sparse, inner format:
 *        \f$ a_{0} + a_{1}*x^{1} + \ldots + a_{n}*x^{n-1} \f$
 */
template<typename T>
class polynomial
{
    static constexpr double pi_23 = (2. * ADF_PI) / 3.;
    const double srt_32 = std::sqrt(3.) / 2.;

    using vec = std::vector<T>;
    using vec_ref = std::vector<T>&;
    using vec_comp = std::vector<complex<T>>;
    using vec_comp_ref = std::vector<complex<T>>&;
    using pair_vec = std::pair<vec, vec>;

    vec m_data;
    bool m_reversed = false;

    /*!
     * \brief _normalize
     * \param other
     */
    void _normalize(vec_ref other)
    {
        for(auto iter = other.rbegin(); iter != other.rend(); ++iter) {
             if(*iter == T(0)) {
                 other.pop_back();
             } else {
                 break;
             }
         }
    }

    /*!
     * \brief _reverse
     * \param other
     */
    void _reverse(vec_ref other)
    {
        std::reverse(other.begin(), other.end());
    }

    /*!
     * \brief _padding
     * \param lhs
     * \param rhs
     */
    void _padding(polynomial<T>& lhs, polynomial<T>& rhs)
    {
        auto add_zero = [&](std::size_t len, std::vector<T> v) -> void
        {
            while(v.size() < len){
                v.push_back(T(0));
            }
        };

        auto lhs_size = lhs.size();
        auto rhs_size = rhs.size();

        if((lhs_size == rhs_size) && lhs_size == 1) {return;}

        else if((lhs_size == rhs_size) && (lhs_size % 2 == 0)) {return;}

        else{
            auto max_size = std::max(lhs_size, rhs_size);

            if(max_size == 2){
                if(lhs_size > rhs_size)
                    add_zero(2, rhs.m_data);
                else
                    add_zero(2, lhs.m_data);
            }
            else if(max_size == 3 || max_size == 4){
                add_zero(4, lhs.m_data);
                add_zero(4, rhs.m_data);
            }
            else{
                add_zero(max_size, lhs.m_data);
                add_zero(max_size, rhs.m_data);
            }
        }
    }

    /*!
     * \brief _shift_pow10
     * \param n
     * \param v
     */
    void _shift_pow10(std::size_t n, vec_ref v)
    {
        std::size_t i = 0;
        while(i < n){
            v.insert(v.begin(), 0);
            ++i;
        }
    }

    /**
     * @brief _pld polynomial long division n / d == q && n % d == r
     * @param n numerator
     * @param d denominator
     * @return std::pair ret - with two vectors = ret.first - quiet
     *         ret.second - remainder
     */
    pair_vec _pld(vec_ref n, vec_ref d)
    {
        if(n.size() == 0 && d.size() == 0)
        {
            return {n, d};
        }

        auto d_size = d.size();
        auto n_degree = n.size() - 1;
        auto d_degree = d.size() - 1;

        auto _q_degree = n_degree - d_degree;
        auto _r_degree = n_degree - d_degree;
        std::size_t _d_tmp;

        vec _d, _q, _r;
        _d.resize(n_degree + 1);
        _q.resize(_q_degree + 1);
        _r.resize(_r_degree + 1);

        if(n_degree >= d_degree)
        {
            while(n_degree >= d_degree)
            {
                _d.assign(d_size, 0);

                for(auto i = 0; i <= d_degree; ++i)
                {
                    _d[i + n_degree - d_degree] = d[i];
                }

                _d_tmp = n_degree;
                _q[n_degree - d_degree] = n[n_degree] / _d[_d_tmp];

                for(auto j = 0; j < _q_degree + 1; ++j)
                {
                    _d[j] = _d[j] * _q[n_degree - d_degree];
                }

                for(auto k = 0; k < n_degree; ++k)
                {
                    n[k] -= _d[k];
                }
                --n_degree;
            }
        }
        for(auto i = 0; i <= n_degree; ++i)
        {
            _r[i] = n[i];
        }
        return {_q, _r};
    }

    /*!
     * \brief polynomial
     * \param v
     */
    polynomial(vec v)
        :m_data(std::move(v))
    {
        _normalize(m_data);
    }

public:
    /**
     * @class polynomial default ctor with initialize order = 0, and not allocate data
     */
    polynomial()
    {}
    /**
     * @brief polynomial for data using input vector
     * @param order first initialize value
     * @param data input array
     *        \f$ a_{n+1}*x^{n} + a_{n}*x^{n-1} + \ldots + a_{1}*x^{1} + a_{0} \f$
     */
    explicit polynomial(const std::initializer_list<T>& list)
        :m_data(std::move(list))
    {
        this->_reverse(m_data);
        m_reversed = true;
    }

    /**
     * @brief polynomial - reserve data for order size
     * @param order - order value
     */
    polynomial(std::size_t size) noexcept
    {
        m_data.resize(size);
        std::fill(m_data.begin(), m_data.end(), 0);
        m_reversed = true;
    }

    /**
     * @brief polynomial - monom implementation
     * @param data - value for fill data array
     */
    explicit polynomial(const T& data) noexcept
    {
        m_data.reserve(1);
        m_data.push_back(std::move(data));
    }

    /**
     * @brief polynomial copy ctor
     * @param other constant object polynomial type
     */
    polynomial(const polynomial<T>& other)
    {
        if(this != &other){
            m_data.reserve(other.size());
            std::copy(other.m_data.begin(), other.m_data.end(), m_data.begin());
        }
    }

    /**
     * @brief operator = copy assignment operator
     * @param other constant object polynomial type
     * @return
     */
    polynomial<T> &operator=(const polynomial<T>& other)
    {
        if(this != &other){
            m_data.clear();
            m_data.reserve(other.size());
            std::copy(other.begin(), other.end(), m_data.begin());
        }
        return (*this);
    }

    /**
     * @brief polynomial move ctor
     * @param other constant object polynomial rvalue type
     */
    polynomial(const polynomial<T>&& other)
    {
        if(this != &other){
            m_data = std::move(other.m_data);
        }
    }

    /**
     * @brief operator = move assignment operator
     * @param other constant object polynomial rvalue type
     * @return
     */
    polynomial<T> &operator=(const polynomial<T>&& other)
    {
        if(this != &other){
            m_data.clear();
            m_data.reserve(other.size());
            m_data = std::move(other.m_data);
        }
        return (*this);
    }

    /**
     * @brief operator = assignment operator
     * @param data initialize polynomial with with one data element
     * @return
     */
    polynomial<T> &operator=(const T& data)
    {
        return *this = polynomial<T>(data);
    }
    
    /*!
     * \brief normalize
     */
    void normalize() { _normalize(m_data); }

    /*!
     * \brief size
     * \return
     */
    std::size_t size() {return m_data.size();}

    /*!
     * \brief size
     * \return
     */
    std::size_t size() const {return m_data.size();}

    /*!
     * \brief empty
     * \return
     */
    bool empty() { return size() == 0;}

    /*!
     * \brief order
     * \return
     */
    std::size_t order()
    {
        return m_data.size() - 1;
    }

    /*!
     * \brief reverse
     * \return
     */
    polynomial<T> reverse()
    {
        polynomial<T> tmp(*this);
        _reverse(tmp.data);
        m_reversed = true;
        return tmp;
    }

    /*!
     * \brief unreverse
     * \return
     */
    polynomial<T> unreverse()
    {
        polynomial<T> tmp(*this);
        _reverse(tmp.m_data);
        m_reversed = false;
        return tmp;
    }

    /*!
     * \brief data
     * \return
     */
    vec data() const
    {
        polynomial<T> tmp = unreverse();
        return tmp.m_data;
    }

    /*!
     * \brief padding
     * \param other
     */
    void padding(polynomial<T>& other) { _padding(*this, other); }

    /*!
     * \brief shift_pow10
     * \param n
     */
    void shift_pow10(std::size_t n) { _shift_pow10(n, this->m_data); }

    /*!
     * \brief roots
     * \param re
     * \param im
     * \return
     */
    int roots(vec_ref re, vec_ref im) const;

    /*!
     * \brief roots
     * \param roots
     * \return
     */
    int roots(vec_comp_ref roots) const;

    /**
     * @brief overloaded operators
     */
    polynomial<T> operator+() { return *this; }
    polynomial<T> operator+() const { return *this; }
    polynomial<T> operator-() const;

    polynomial<T> operator+(const polynomial<T>& other) const;
    polynomial<T> operator-(const polynomial<T>& other) const;
    polynomial<T> operator*(const polynomial<T>& other) const;
    polynomial<T> operator/(const polynomial<T>& other) const;
    polynomial<T> operator%(const polynomial<T>& other) const;

    polynomial<T> operator^(uint32_t n) const;

    polynomial<T> operator+(const T& data) const;
    polynomial<T> operator-(const T& data) const;
    polynomial<T> operator*(const T& data) const;
    polynomial<T> operator/(const T& data) const;
    polynomial<T> operator%(const T& data) const;

    polynomial<T>& operator+=(const polynomial<T>& other);
    polynomial<T>& operator-=(const polynomial<T>& other);
    polynomial<T>& operator*=(const polynomial<T>& other);
    polynomial<T>& operator/=(const polynomial<T>& other);
    polynomial<T>& operator%=(const polynomial<T>& other);

    polynomial<T>& operator+=(const T& data);
    polynomial<T>& operator-=(const T& data);
    polynomial<T>& operator*=(const T& data);
    polynomial<T>& operator/=(const T& data);
    polynomial<T>& operator%=(const T& data);

    T& operator[](std::size_t i) {return this->m_data.at(i);}
    const T& operator[](std::size_t i) const {return this->m_data.at(i);}

    /**
     * @brief out stream
     */
    template<typename U>
    friend std::ostream &operator<<(std::ostream &p_out, const polynomial<U> &p_val);

private:
    /*!
     *
     */
    template<typename U>
    friend polynomial<U> karatsuba(polynomial<U>& rhs,
                                   polynomial<U>& lhs);

    //TODO
    /*!
     *
     */
    template<typename U>
    friend polynomial<U> newton(polynomial<U>& rhs,
                                polynomial<U>& lhs);

    /*!
     * \brief quadratic
     * \param a
     * \param b
     * \param c
     * \param re1
     * \param im1
     * \param re2
     * \param im2
     * \return
     */
    int quadratic(const T& a, const T& b, const T& c,
                  T& re1, T& im1,
                  T& re2, T& im2) const;

    /*!
     * \brief quadratic
     * \param a
     * \param b
     * \param c
     * \param x1
     * \param x2
     * \return
     */
    int quadratic(const T& a, const T& b, const T& c,
                  complex<T>& x1, complex<T>& x2) const;

    /*!
     * \brief quadratic
     * \param b
     * \param c
     * \param re1
     * \param im1
     * \param re2
     * \param im2
     * \return
     */
    int quadratic(const T& b, const T& c,
                  T& re1, T& im1,
                  T& re2, T& im2) const;

    /*!
     * \brief quadratic
     * \param b
     * \param c
     * \param x1
     * \param x2
     * \return
     */
    int quadratic(const T& b, const T& c,
                  complex<T>& x1, complex<T>& x2) const;

    /*!
     * \brief cubic
     * \param a
     * \param b
     * \param c
     * \param re1
     * \param im1
     * \param re2
     * \param im2
     * \param re3
     * \param im3
     * \return
     */
    int cubic(const T& a, const T& b, const T& c,
              T& re1, T& im1,
              T& re2, T& im2,
              T& re3, T& im3) const;

    /*!
     * \brief cubic
     * \param a
     * \param b
     * \param c
     * \param x1
     * \param x2
     * \param x3
     * \return
     */
    int cubic(const T& a, const T& b, const T& c,
              complex<T>& x1,
              complex<T>& x2,
              complex<T>& x3) const;

    /*!
     * \brief laguerre
     * \param poly
     * \param degree
     * \param root
     * \return
     */
    int laguerre(vec_comp_ref poly,
                 const std::size_t degree,
                         complex<T>& root) const;

    /*!
     * \brief find_roots
     * \param poly
     * \param roots
     * \param polish
     * \return
     */
    int find_roots(vec_ref poly,
                   vec_comp_ref roots,
                   const bool polish = true) const;

};

template<typename U>
std::ostream &operator<<(std::ostream &p_out, const polynomial<U> &p_val)
{
    if(p_val.empty()) return p_out << U(0);
    polynomial<U> tmp(p_val);
    tmp.unreverse();
    for(auto i=0; i<tmp.size(); ++i)
    {
        p_out << tmp[i];
    }
    return p_out;
}

template<typename U>
polynomial<U> karatsuba(polynomial<U>& lhs, polynomial<U>& rhs)
{
    auto lhs_size = lhs.size();
    auto rhs_size = rhs.size();

    if(lhs_size == 1 && rhs_size == 1)
    {
        polynomial<U> ans(lhs_size);
        U m = lhs[0] * rhs[0];
        ans.m_data[0] = m;
        return ans;
    }
    else{

        lhs.padding(rhs);
        auto len = lhs.size();
        std::size_t half = len / 2;

        polynomial<U> result(len * 2);

        /* Split the digit sequences in the middle. */
        polynomial<U> a_lhs(half);
        polynomial<U> a_rhs(half);

        for(std::size_t i = 0; i < half; ++i)
        {
            a_lhs.m_data[i] = lhs.m_data[i];
            a_rhs.m_data[i] = rhs.m_data[i];
        }

        polynomial<U> b_lhs(half);
        polynomial<U> b_rhs(half);

        for(std::size_t i = half; i < len; ++i)
        {
            b_lhs.m_data[i - half] = lhs.m_data[i];
            b_rhs.m_data[i - half] = rhs.m_data[i];
        }

        /* 3 recursive calls made to numbers approximately half the size. */
        polynomial<U> c_2 = karatsuba(b_lhs, b_rhs);
        polynomial<U> c_0 = karatsuba(a_lhs, a_rhs);

        polynomial<U> a_mid = b_lhs + a_lhs;
        polynomial<U> b_mid = b_rhs + a_rhs;

        polynomial<U> p_1 = karatsuba(a_mid, b_mid);
        polynomial<U> tmp = c_2 + c_0;
        polynomial<U> c_1 = p_1 - tmp;

        if((len % 2) != 0){
            c_2.shift_pow10(len - 1);
        }
        else{
            c_2.shift_pow10(len);
        }
        c_1.shift_pow10(half);

        polynomial<U> c_21 = c_2 + c_1;
        c_21.normalize();
        result = c_21 + c_0;

        return result;
    }
}

template<typename U>
polynomial<U> newton(const polynomial<U>& rhs,
                     const polynomial<U>& lhs)
{
    U two = one(U()) + one(U());
    polynomial<U> g = one(U());
    polynomial<U> f = lhs;

    auto n = (rhs.empty()) ? 0 : rhs.order();
    auto m = (lhs.empty()) ? 0 : lhs.order();

    int i, j, r;

    if(n < m) return polynomial<U>(zero(U()));

    U lc = (lhs.empty()) ? one(U()) : lhs.m_data.front();
    U lci = one(U()) / lc;
    f *= lci;

    n = n - m + 1;
    r = 0;
    j = 1;
    while(j < n) {r++; j *= 2;}

    for(i = 1, j = 2; i <= r; i++, j *= 2)
    {
        g *= (two - g * f);
        while(!g.empty() && g.order() >= j)
            g.m_data.pop_front();
    }

    polynomial<U> result = (rhs.reverse()) * g * lci;
    while(!result.empty() && result.order() >= n)
        result.m_data.pop_front();
    return result.reverse();
}

template<typename T>
int polynomial<T>::quadratic(const T& a, const T& b, const T& c,
                             T& re1, T& im1,
                             T& re2, T& im2) const
{
    const T aa = (a + a);
    const T det = (b * b) - (2 * aa * c);

    if(det >= 0)
    {
        im1 = im2 = 0.;
        const T sq_det = std::sqrt(det);
        re1 = (-b + sq_det) / aa;
        re2 = (-b - sq_det) / aa;
        return ((det == 0.0f) ? 1 : 2);
    }
    else{
        const T sq_det = std::sqrt(-det);
        re1 = re2 = -b / aa;
        im1 = sq_det;
        im2 = -sq_det;
    }
    return 0;
}

template<typename T>
int polynomial<T>::quadratic(const T& a, const T& b, const T& c,
                             complex<T>& x1, complex<T>& x2) const
{
    const complex<T> det = b * b - 4 * a * c;
    const complex<T> q = -0.5 * (b + sgn(b) * std::sqrt(det));
    x1 = q/a;
    x2 = c/q;

    return (det < 0.0) ? 0 : ((det == 0.0) ? 1 : 2);
}

template<typename T>
int polynomial<T>::quadratic(const T& b, const T& c,
                             T& re1, T& im1,
                             T& re2, T& im2) const
{
    const T hb = (b / 2.);
    const T det = hb * hb - c;

    if(det >= 0)
    {
        im1 = im2 = 0.0;
        const T sq_det = std::sqrt(det);
        re1 = -hb + sq_det;
        re2 = -hb - sq_det;
        return ((det == 0.0f) ? 1 : 2);
    }
    else {
        const T sq_det = std::sqrt(-det);
        re1 = re2 = -hb;
        im1 = sq_det;
        im2 = -sq_det;
    }
    return 0;
}

template<typename T>
int polynomial<T>::quadratic(const T& b, const T& c,
                             complex<T>& x1, complex<T>& x2) const
{
    const T det = b * b - 4 * c;
    const complex<T> cc = -0.5 * (b + sgn(b) * std::sqrt(det));
    x1 = cc;
    x2 = c / cc;

    return (det < 0.0) ? 0 : ((det == 0.0) ? 1 : 2);
}

template<typename T>
int polynomial<T>::cubic(const T& a, const T& b, const T& c,
                         T& re1, T& im1,
                         T& re2, T& im2,
                         T& re3, T& im3) const
{
    complex<T> x1, x2, x3;
    int out = cubic(a, b, c, x1, x2, x3);

    re1 = x1.getReal();  im1 = x1.getImag();
    re2 = x2.getReal();  im2 = x2.getImag();
    re3 = x3.getReal();  im3 = x3.getImag();

    return out;
}

template<typename T>
int polynomial<T>::cubic(const T& a, const T& b, const T& c,
                         complex<T>& x1,
                         complex<T>& x2,
                         complex<T>& x3) const
{
    T tmp = (a * a);
    const T q1 = ((tmp - 3.0 * b) / 9.0);
    const T r1 = ((2.0 * tmp * a - 9.0 * a * b + 27.0 * c) / 54.0);
    const T r2 = r1 * r1;
    const T q2 = q1 * q1 * q1;
    tmp = a / 3.0;

    if(r2 < q2)
    {
        T q_sqrt = std::sqrt(q1);
        const T theta = std::acos(r1 / (q1 * q_sqrt)) / 3.0;
        q_sqrt *= -2.0;
        x1 = q_sqrt * std::cos(theta ) - tmp;
        x2 = q_sqrt * std::cos(theta + pi_23) - tmp;
        x3 = q_sqrt * std::cos(theta - pi_23) - tmp;
        return 3;
    }
    else {
        constexpr double aa = -std::pow(r1 + sgn(r1) * std::sqrt(r2 - q2), 1.0 / 3.0);
        T bb = 0.0;
        if(aa != 0.0)
        {
            bb = q1 / aa;
        }
        const T apb = aa + bb;
        const T amb = aa - bb;

        x1 = apb - tmp;
        x2 = complex<T>(-0.5 * apb - tmp, srt_32 * amb);
        x3 = complex<T>(-0.5 * apb - tmp, -srt_32 * amb);

        return 1;
    }
}

template<typename T>
int polynomial<T>::laguerre(vec_comp_ref poly,
                            const std::size_t degree,
                            complex<T>& root) const
{
    static const double epsilon = std::numeric_limits<double>::epsilon();

    T absx, absp, absm, err;
    complex<T> dx, x1, b, d, f, g, h, sq, gp, gm, g2;

    static const double frac[] = {0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 1.};

    const int mm1 = degree - 1;
    int it;

    for(it = 1; it <= MAX_IT; ++it)
    {
        b = poly.at(degree);
        err = abs(b);
        absx = abs(root);

        for(auto j = mm1; j >=0; --j)
        {
            f = root * f + d;
            d = root * d + b;
            b = root * b + poly.at(j);
            err = b.abs() + absx * err;
        }

        err *= epsilon;

        if(abs(b) < err) return it;

        //Laguerre method
        g = b / b; // G(x_k) = P'(x_k) / P(x_k). ([1], eq. 7.10)
        g2 = g * g;
        h = g2 - 2. * f / b; // H(x_k) = G(x_k)^2 - (P''(x_k)/P(x_k)). ([1], eq. 7.11)

        auto tmp = static_cast<T>(mm1) * complex<T>(static_cast<T>(degree)) * h - g2;
        sq = sqrt(tmp);
        gp = g + sq;
        gm = g - sq;
        absp = abs(gp);
        absm = abs(gm);

        if(absp < absm) gp = gm;

        auto cdegree = complex<T>(static_cast<T>(degree)) / gp;
        auto cabsx = (1. + absx) * complex<T>(std::cos(T(it)), std::sin(T(it)));
        dx = ((std::max(absp, absm) > 0.0) ? cdegree : cabsx);
        x1 = root - dx;

        if(root == x1) return it;

        if((it % MT) != 0)
        {
            root = x1;
        }
        else {
            root = root - frac[it / MT] * dx;
        }
    }
    return it;
}

template<typename T>
int polynomial<T>::find_roots(vec_ref poly,
                              vec_comp_ref roots,
                              const bool polish) const
{
    auto vec_fill = [&](vec_ref other) -> vec_comp
    {
        vec_comp out;
        out.resize(other.size(), complex<T>(0., 0.));

        typename vec::const_iterator it, eit;
        typename vec_comp::iterator dit;

        for(it = other.begin(), dit = out.begin(), eit = other.end();
            it!=eit; ++it, ++eit)
        {
            (*dit).setReal(static_cast<T>(*it));
        }
        return out;
    };

    static const double epsilon = std::numeric_limits<double>::epsilon();

    const int degree = poly.size() - 1;
    roots.resize(degree);

    vec_comp ad = vec_fill(poly);

    int i, its, j, jj;
    complex<T> x, b, c;

    for(j = degree; j >= 1; --j)
    {
        laguerre(ad, j, x);

        if(std::abs(x.getImag()) <= 2. * epsilon * std::abs(x.getReal()))
        {
            x = complex<T>(x.getReal(), 0);
        }

        roots.at(j - 1) = x;

        b = ad.at(j);
        for(jj = j-1; jj >= 0; --jj)
        {
            c = ad.at(jj);
            ad.at(jj) = b;
            b = (x * b) + c;
        }
    }

    if(polish)
    {
        ad = vec_fill(poly);
        for(j = 0; j < degree; ++j)
        {
            laguerre(ad, degree, roots[j]);
        }
    }

    its = (roots.at(0).getImag() == 0.0) ? 1 : 0;
    for(j = 1; j < degree; ++j)
    {
        x = roots.at(j);
        if(x.getImag() == 0.0) its++;

        for(i = j - 1; i >= 0; --i)
        {
            if(roots.at(i).getReal() <= x.getReal()) break;

            roots.at(i + 1) = roots.at(i);
        }
        roots.at(i + 1) = x;
    }
    return its;
}

template<typename T>
int polynomial<T>::roots(vec_ref re, vec_ref im) const
{
    const std::size_t indx = m_data.size() - 1;
    re.resize(indx, zero(T()));
    im.resize(indx, zero(T()));

    switch(indx)
    {
    case 1:
        re.at(0) = -m_data.at(0)/m_data.at(1);
        im.at(0) = 0.;
        return 1;
    case 2:
        return quadratic(m_data.at(2), m_data.at(1), m_data.at(0),
                         re.at(0), im.at(0),
                         re.at(1), im.at(1));
    case 3:
        if(m_data.at(3) == zero(T()))
        {
            return quadratic(m_data.at(2), m_data.at(1), m_data.at(0),
                             re.at(0), im.at(0),
                             re.at(1), im.at(1));
        } else {
            return cubic(m_data.at(2)/m_data.at(3), m_data.at(1)/m_data.at(3), m_data.at(0)/m_data.at(3),
                         re.at(0), im.at(0),
                         re.at(1), im.at(1),
                         re.at(2), im.at(2));
        }
        break;
    default:
        vec_comp roots_out;
        int out = find_roots(m_data, roots_out);
        for(auto i=0; i<roots_out.size(); ++i)
        {
            re.at(i) = roots_out.at(i).getReal();
            im.at(i) = roots_out.at(i).getImag();
        }
        return out;
    }
    return -1;
}

template<typename T>
int polynomial<T>::roots(vec_comp_ref roots) const
{
    const std::size_t indx = m_data.size() - 1;
    roots.resize(indx, complex<T>(0., 0.));

    switch(indx)
    {
    case 1:
        roots.at(0) = complex<T>(-m_data.at(0) / m_data.at(1));
        return 1;
    case 2:
        return quadratic(m_data.at(2), m_data.at(1), m_data.at(0),
                         roots.at(0), roots.at(1));
    case 3:
        if(m_data.at(3) == 0.0)
        {
            return quadratic(m_data.at(2), m_data.at(1), m_data.at(0),
                             roots.at(0), roots.at(1));
        } else {
            return cubic(m_data.at(2) / m_data.at(3), m_data.at(1) / m_data.at(3), m_data.at(0) / m_data.at(3),
                         roots.at(0), roots.at(1), roots.at(2));
        }
        break;
    default:
        return find_roots(m_data, roots);
    }
    return -1;
}

template<typename U> bool operator==(const polynomial<U> &p_fr, const polynomial<U> &p_sc)
{
    if(p_fr.size() != p_sc.size())
        return false;
    for(std::size_t i = 0; i < p_fr.size(); ++i)
    {
        if(p_fr[i] != p_sc[i])
            return false;
    }
    return true;
}

template<typename U> bool operator!=(const polynomial<U> &p_fr, const polynomial<U> &p_sc)
{
    return p_fr == p_sc;
}

/**
 * @brief implementation
 */
template<typename T>
polynomial<T> polynomial<T>::operator-() const
{
    polynomial<T> tmp(*this);
    for(auto i=0; i<tmp.size(); ++i) tmp[i] = -(tmp[i]);
    return tmp;
}

template<typename T>
polynomial<T> polynomial<T>::operator+(const polynomial<T>& other) const
{
    std::size_t ms = std::max(this->m_data.size(), other.size());
    polynomial<T> temp(ms);

    for(int i = ms-1; i >= std::min(this->m_data.size(), other.size()); --i)
    {
        temp[i] = (this->m_data.size(), other.size() ? this->m_data[i] : other[i]);
    }
    for(int j = std::min(this->m_data.size(), other.size())-1; j >=0; --j)
    {
        temp[j] = this->m_data[j] + other[j];
    }
    temp.normalize();
    return temp;
}

template<typename T>
polynomial<T> polynomial<T>::operator-(const polynomial<T>& other) const
{
    return (*this) + (-other);
}

template<typename T>
polynomial<T> polynomial<T>::operator*(const polynomial<T> &other) const
{
    polynomial<T> result = karatsuba(*this, other);
    result.normalize();
    return result;
}

template<typename T>
polynomial<T> polynomial<T>::operator/(const polynomial<T> &other) const
{
    pair_vec result = _pld(*this->m_data, other.m_data);
    return polynomial(result.first);
}

template<typename T>
polynomial<T> polynomial<T>::operator%(const polynomial<T>& other) const
{
    pair_vec result = _pld(*this->m_data, other.m_data);
    return polynomial(result.second);
}

template<typename T>
polynomial<T> polynomial<T>::operator^(uint32_t n) const
{
    polynomial<T> result(one<T>(T())), factor(*this);
    while(n > 0)
    {
        if(n % 2 == 1) result *= factor;
        factor *= factor;
        n /= 2;
    }
    return result;
}

template<typename T>
polynomial<T> polynomial<T>::operator+(const T& data) const
{
    return (*this) + polynomial<T>(data);
}

template<typename T>
polynomial<T> polynomial<T>::operator-(const T& data) const
{
    return (*this) - polynomial<T>(data);
}

template<typename T>
polynomial<T> polynomial<T>::operator*(const T& data) const
{
    return (*this) - polynomial<T>(data);
}

template<typename T>
polynomial<T> polynomial<T>::operator/(const T& data) const
{
    return (*this) / polynomial<T>(data);
}

template<typename T>
polynomial<T> polynomial<T>::operator%(const T& data) const
{
    return (*this) / polynomial<T>(data);
}

template<typename T>
polynomial<T>& polynomial<T>::operator+=(const polynomial<T>& other)
{
    return *this = *this + other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator-=(const polynomial<T>& other)
{
    return *this = *this - other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator*=(const polynomial<T>& other)
{
    return *this = *this * other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator/=(const polynomial<T>& other)
{
    return *this = *this / other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator%=(const polynomial<T>& other)
{
    return *this = *this % other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator+=(const T& data)
{
    return *this += polynomial<T>(data);
}

template<typename T>
polynomial<T>& polynomial<T>::operator-=(const T& data)
{
    return *this -= polynomial<T>(data);
}

template<typename T>
polynomial<T>& polynomial<T>::operator*=(const T& data)
{
    return *this *= polynomial<T>(data);
}

template<typename T>
polynomial<T>& polynomial<T>::operator/=(const T& data)
{
    return *this /= polynomial<T>(data);
}

template<typename T>
polynomial<T>& polynomial<T>::operator%=(const T& data)
{
    return *this %= polynomial(data);
}

} //adf

#endif //COMPLEX_H
