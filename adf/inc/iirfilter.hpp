// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

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
