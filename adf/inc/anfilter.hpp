#ifndef ANALOGFILTER_H
#define ANALOGFILTER_H

#include "genfilter.hpp"

namespace adf {

template<typename T=double>
class AnalogFilter
{
    CalcCoeffs<T>* m_calccoeffs;

public:
    explicit AnalogFilter<T>(FiltParam<T>& fparam, FilterType& ftype, ApproxType& type)
    {
        m_calccoeffs = new CalcCoeffs<T>(fparam, ftype, type);
    }

    ~AnalogFilter<T>()
    {
        delete m_calccoeffs;
    }


};

} //adf

#endif //ANALOGFILTER_H
