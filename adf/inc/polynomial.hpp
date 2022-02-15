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
    std::size_t m_order;

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

public:
    /**
     * @class polynomial default ctor with initialize order = 0, and not allocate data
     */
    polynomial() noexcept
    {
        m_order=0;
    }
    /**
     * @brief polynomial for data using input vector
     * @param order first initialize value
     * @param data input array
     *        \f$ a_{n+1}*x^{n} + a_{n}*x^{n-1} + \ldots + a_{1}*x^{1} + a_{0} \f$
     */
    explicit polynomial(std::vector<T> data, std::size_t order) noexcept
    {
        assert(((m_order == 0)||(!data.empty())) && " Zero order or array empty");
        assert((m_order == (data.size()+1)) && "Order and data size+1 not equal");
        m_order = order;
        m_data.reserve(m_order+1);
        m_data = std::move(data);
    }

    /**
     * @brief polynomial - reserve data for order size
     * @param order - order value
     */
    explicit polynomial(std::size_t order) noexcept
        :m_order(order)
    {
        m_data(order, 0);
    }

    /**
     * @brief polynomial - monom implementation
     * @param data - value for fill data array
     */
    explicit polynomial(const T& data) noexcept
    {
        m_order = 1;
        m_data.reserve(m_order);
        m_data.push_back(data);
    }

    /**
     * @brief polynomial copy ctor
     * @param other constant object polynomial type
     */
    polynomial(const polynomial<T>& other)
        :m_order(other.m_order)
    {
        m_data.reserve((other.m_order)+1);
        std::copy(other.begin(), other.end(), m_data.begin());
    }

    /**
     * @brief operator = copy assignment operator
     * @param other constant object polynomial type
     * @return
     */
    polynomial<T> &operator=(const polynomial<T>& other)
    {
        if(this != &other){
            m_order = other.m_order;
            m_data.reserve((other.m_order)+1);
            std::copy(other.begin(), other.end(), m_data.begin());
        }
        return (*this);
    }

    /**
     * @brief polynomial move ctor
     * @param other constant object polynomial rvalue type
     */
    polynomial(const polynomial<T>&& other)
        :m_order(std::move(other.m_data))
    {
        m_data.reserve((other.m_order)+1);
        m_data = std::move(other.m_data);
    }

    /**
     * @brief operator = move assignment operator
     * @param other constant object polynomial rvalue type
     * @return
     */
    polynomial<T> &operator=(const polynomial<T>&& other)
    {
        if(this != &other){
            m_order = std::move(other.m_order);
            m_data.reserve((other.m_order)+1);
            m_data = std::move(other.m_data);
        }
        return (*this);
    }
    
    explicit polynomial(const std::initializer_list<T>& list)
        :m_data(list)
        ,m_order(list.size())
    {}

    auto normalize() -> void { _normalize(m_data); }
    auto order() -> std::size_t { return this->m_order; }
    auto reverse() -> polynomial<T>;

    /**
     * @brief overloaded operators
     */
    polynomial<T> operator+() { return *this; }
    polynomial<T> operator+() const { return *this; }
    polynomial<T> operator-() const;

    polynomial<T> operator+(polynomial<T>& other);
    polynomial<T> operator+(polynomial<T>& other) const;
    polynomial<T> operator-(polynomial<T>& other);
    polynomial<T> operator-(polynomial<T>& other) const;
    polynomial<T> operator*(polynomial<T>& other);
    polynomial<T> operator*(polynomial<T>& other) const;
    polynomial<T> operator/(polynomial<T>& other);
    polynomial<T> operator/(polynomial<T>& other) const;
    polynomial<T> operator%(polynomial<T>& other);
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
    T& operator[](std::size_t i) const {return this->m_data[i];}

    /**
     * @brief out stream
     */
    template<typename U>
    friend std::ostream &operator<<(std::ostream &p_out, const polynomial<U> &p_val);
};

template<typename U> polynomial<U> operator+(polynomial<U> &p_fr, polynomial<U> &p_sc)
{
    auto ms = std::max(p_fr.order(), p_sc.order());
    polynomial<U> temp(ms);

    for(auto i = ms-1; i >= std::min(p_fr.order(), p_sc.order()); --i)
    {
        temp[i] = (p_fr.order() > p_sc.order() ? p_fr[i] : p_sc[i]);
    }
    for(auto j = std::min(p_fr.order(), p_sc.order())-1; j >=0; --j)
    {
        temp[j] = p_fr[j] + p_sc[j];
    }
    temp.normalize();
    return temp;
}

template<typename U> polynomial<U> operator-(polynomial<U> &p_fr, polynomial<U> &p_sc)
{

}

template<typename U> polynomial<U> operator*(polynomial<U> &p_fr, polynomial<U> &p_sc)
{

}

template<typename U> polynomial<U> operator/(polynomial<U> &p_fr, polynomial<U> &p_sc)
{

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
    for(auto i=0; i<tmp.m_data.size(); ++i) tmp.m_data[i] = -(tmp.m_data);
    return tmp;
}

template<typename T>
polynomial<T> polynomial<T>::operator+(polynomial<T>& other)
{
    auto ms = std::max(this->order(), other.order());
    polynomial<T> temp(ms);

    for(auto i = ms-1; i >= std::min(this->order(), other.order()); --i)
    {
        temp[i] = (this->order() > other.order() ? this->m_data[i] : other[i]);
    }
    for(auto j = std::min(this->order(), other.order())-1; j >=0; --j)
    {
        temp[j] = this->m_data[j] + other[j];
    }
    temp.normalize();
    return temp;
}

template<typename T>
polynomial<T> polynomial<T>::operator+(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator-(polynomial<T>& other)
{

}

template<typename T>
polynomial<T> polynomial<T>::operator-(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator*(polynomial<T>& other)
{

}

template<typename T>
polynomial<T> polynomial<T>::operator*(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator/(polynomial<T>& other)
{

}

template<typename T>
polynomial<T> polynomial<T>::operator/(polynomial<T>& other) const
{

}

template<typename T>
polynomial<T> polynomial<T>::operator%(polynomial<T>& other)
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
