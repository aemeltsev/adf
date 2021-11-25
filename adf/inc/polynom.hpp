#ifndef COMPLEX_H
#define COMPLEX_H

#include <vector>
#include <algorithm>

#include <cassert>

namespace adf
{
/**
 * @class Polynomial - parameterized class for polynomial working
 *        inner format:
 *        \f$ a_{n+1}*x^{n} + a_{n}*x^{n-1} + \ldots + a_{1}*x^{1} + a_{0} \f$
 */
template<typename T>
class polynomial
{
    std::vector<T> m_data;
    std::size_t m_order;

    void normalize(std::vector<T>& other)
    {
        for(auto iter = other.rbegin(); iter != other.rend(); ++iter) {
             if (*iter == T(0)) {
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
     * @class polynomial for data using input vector
     * @param order first initialize value
     * @param data input array
     *        \f$ a_{n+1}*x^{n} + a_{n}*x^{n-1} + \ldots + a_{1}*x^{1} + a_{0} \f$
     */
    polynomial(std::size_t order, std::vector<T>& data) noexcept
        :m_order(order)
    {
        assert(((m_order == 0)||(!data.empty())) && " Zero order or array empty");
        assert((m_order == (data.size()+1)) && "Order and data size+1 not equal");
        m_data.reserve(m_order+1);
        m_data = std::move(data);
    }

    /**
     * @class polynomial copy ctor
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
    
};

} //adf

#endif //COMPLEX_H
