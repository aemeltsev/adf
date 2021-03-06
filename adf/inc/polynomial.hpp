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

namespace adf
{

/**
 * @class Polynomial - parameterized class for polynomial working
 *        not sparse, inner format:
 *        \f$ a_{0} + a_{1}*x^{1} + \ldots + a_{n}*x^{n-1} \f$
 */
template<typename T>
class polynomial
{
    using vec = std::vector<T>;
    using vec_ref = std::vector<T>&;
    using pair_vec = std::pair<vec, vec>;

    vec m_data;
    bool m_reversed = false;

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

    void _reverse(vec_ref other)
    {
        std::reverse(other.begin(), other.end());
    }

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
    
    void normalize() { _normalize(m_data); }

    std::size_t size() {return m_data.size();}

    std::size_t size() const {return m_data.size();}

    bool empty() { return size() == 0;}

    polynomial<T> unreverse()
    {
        polynomial<T> tmp(*this);
        _reverse(tmp.m_data);
        m_reversed = false;
        return tmp;
    }

    vec data() const
    {
        polynomial<T> tmp = unreverse();
        return tmp.m_data;
    }

    void padding(polynomial<T>& other) { _padding(*this, other); }

    void shift_pow10(std::size_t n) { _shift_pow10(n, this->m_data); }

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
    template<typename U>
    friend polynomial<U> karatsuba(polynomial<U>& rhs,
                                   polynomial<U>& lhs);

    //TODO
    template<typename U>
    friend polynomial<U> newton(polynomial<U>& rhs,
                                polynomial<U>& lhs);

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
/*
template<typename U>
polynomial<U> newton(const polynomial<U>& rhs,
                     const polynomial<U>& lhs)
{

}
*/
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
