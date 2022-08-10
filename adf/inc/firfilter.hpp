// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef DIGFIRFILTER_H
#define DIGFIRFILTER_H

#include "genfilter.hpp"

namespace adf {

template<typename T>
class DigFIRFilter
{
private:

    FiltParam<T> m_fparam;
    FilterType m_ftype;
    ApproxType m_atype;
    std::size_t m_order = 0;
    std::vector<T> n_acoefs, n_bcoefs; /*!< to normalise coefs */
    std::vector<T> un_acoefs, un_bcoefs; /*!< to unnormalise coefs */
public:

    enum class WindowType
    {
        RectangularWindow,
        HammingWindow,
        HanningWindow,
        BlackmanWindow,
        KaiserWindow,
        ParksMcClellanWindow,
        BartlettWindow
    };

    DigFIRFilter(FilterType& f_type, FiltParam<T>& f_param, ApproxType& atype)
    {

    }

    void CalcDigFIRCoeffs(const WindowType& w_type);
    void CalcFilterLen();
    void CalcIdealFIRCoeffs();
    void RectangularWindowCoeffs();
    void HammingWindowCoeffs();
    void HanningWindowCoeffs();
    void BlackmanWindowCoeffs();
    void KaiserWindowCoeffs();
    void BartlettWindowCoeffs();
    void ParksMcClellanWindowCoeffs();
    void setParksMcClellanParameters();
    void RemezInterchange();
    void ComputeLagrange();
    void EstimateFreqResponse();
    void ComputeCoeffs();
    void MultWinIdealCoeffs();
};

} //namespace adf

#endif //DIGFIRFILTER_H
