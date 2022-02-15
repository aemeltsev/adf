#include "inc/polynomial.hpp"
namespace adf
{

template <typename T>
std::vector<T> conv(const std::vector<double>& lhs, const std::vector<double>& rhs)
{
    auto size = (lhs.size() + rhs.size() - 1);
    std::vector<T> result(size, 0);

    for(std::size_t k=0; k<lhs.size(); ++k){
        for(std::size_t n=0; n<rhs.size(); ++n){
            result[k + n] += lhs[k] * rhs[n];
        }
    }
    return result;
}

}
