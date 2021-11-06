#ifndef DIGIIRFILTER_H
#define DIGIIRFILTER_H

#include "genfilter.hpp"

namespace adf {

template<typename T=double>
class DigIIRFilter
{
    CalcCoeffs<T>* m_calccoeffs;

public:
    explicit DigIIRFilter<T>(FiltParam<T>& fparam, FilterType& ftype, ApproxType& type)
    {
        m_calccoeffs = new CalcCoeffs<T>(fparam, ftype, type);
    }

    ~DigIIRFilter<T>()
    {
        delete m_calccoeffs;
    }


};

}

#endif //DIGIIRFILTER_H
