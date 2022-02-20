/**
  * License text...
*/
#ifndef POLYNOM_H
#define POLYNOM_H

#include <iostream>
#include <cassert>
#include <vector>
#include <utility>
#include <algorithm>
#include <initializer_list>

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
    std::vector<T> m_data;
    bool m_reversed = false;

    void _normalize(std::vector<T>& other)
    {
        for(auto iter = other.rbegin(); iter != other.rend(); ++iter) {
             if(*iter == T(0)) {
                 other.pop_back();
             } else {
                 break;
             }
         }
    }

    void _reverse(std::vector<T>& other)
    {
        std::reverse(other.begin(), other.end());
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
    std::vector<T> data() const
    {
        polynomial<T> tmp = unreverse();
        return tmp.m_data;
    }

    /**
     * @brief overloaded operators
     */
    polynomial<T> operator+() { return *this; }
    polynomial<T> operator+() const { return *this; }
    polynomial<T> operator-() const;

    polynomial<T> operator+(const polynomial<T> &other) const;
    polynomial<T> operator-(polynomial<T>& other) const;
    polynomial<T> operator*(polynomial<T>& other) const;
    polynomial<T> operator/(polynomial<T>& other) const;
    polynomial<T> operator%(polynomial<T>& other) const;

    polynomial<T> operator^(polynomial<T>& other) const;

    polynomial<T> operator+(const T& data) const;
    polynomial<T> operator-(const T& data) const;
    polynomial<T> operator*(const T& data) const;
    polynomial<T> operator/(const T& data) const;
    polynomial<T> operator%(const T& data) const;

    polynomial<T>& operator+=(const polynomial<T>& other);
    polynomial<T>& operator-=(const polynomial<T>& other);
    polynomial<T>& operator*=(const polynomial<T>& other);
    polynomial<T>& operator/=(const polynomial<T>& other);

    polynomial<T>& operator+=(const T& data);
    polynomial<T>& operator-=(const T& data);
    polynomial<T>& operator*=(const T& data);
    polynomial<T>& operator/=(const T& data);

    T& operator[](std::size_t i) {return this->m_data[i];}
    const T& operator[](std::size_t i) const {return this->m_data[i];}

    /**
     * @brief out stream
     */
    template<typename U>
    friend std::ostream &operator<<(std::ostream &p_out, const polynomial<U> &p_val);
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

template<typename U> bool operator==(polynomial<U> &p_fr, polynomial<U> &p_sc);
template<typename U> bool operator!=(polynomial<U> &p_fr, polynomial<U> &p_sc);
template<typename U> bool operator==(const polynomial<U> &p_fr, const polynomial<U> &p_sc);
template<typename U> bool operator!=(const polynomial<U> &p_fr, const polynomial<U> &p_sc);

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
polynomial<T> polynomial<T>::operator-(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator*(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator/(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator%(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator^(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator+(const T& data) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator-(const T& data) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator*(const T& data) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator/(const T& data) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator%(const T& data) const
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator+=(const polynomial<T>& other)
{
    return *this = *this + other;
}

template<typename T>
polynomial<T>& polynomial<T>::operator-=(const polynomial<T>& other)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator*=(const polynomial<T>& other)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator/=(const polynomial<T>& other)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator+=(const T& data)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator-=(const T& data)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator*=(const T& data)
{

}

template<typename T>
polynomial<T>& polynomial<T>::operator/=(const T& data)
{

}

} //adf

#endif //COMPLEX_H