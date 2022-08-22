// Copyright (C) 1994 Les Thede
// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef DIGFIRFILTER_H
#define DIGFIRFILTER_H

#include "advmath.hpp"
#include "genfilter.hpp"

constexpr int ADF_GRID_MULT(16);         /**< starting grid mult */
constexpr int ADF_MAX_BANDS(10);         /**< max bands for PM */
constexpr double ADF_INF(1E+300);
constexpr double ADF_ERR_EPS(1E-6);

namespace adf {

enum class WindowType
{
    RectangularWindow,
    HammingWindow,
    HanningWindow,
    BlackmanWindow,
    KaiserWindow,
    ParksMcClellanWindow,
    BartlettWindow,
    Undefined
};

enum class FIRType
{
    SymOdd = 1,
    SymEven = 2,
    AsymOdd = 3,
    AsymEven = 4,
};

template<typename T>
class DigFIRFilter
{
private:

    FiltParam<T> m_fparam;
    FilterType m_ftype;
    ApproxType m_atype;
    std::size_t m_order = 0;
    std::vector<T> n_acoefs, n_bcoefs;

    struct ParksMcClellanParam
    {
        int num_pts;
        int num_ext;
        int num_bands;
        FIRType filt_type;
        int grid_mult = ADF_GRID_MULT;
        std::vector<int> band_pts;
        std::vector<int> extrml;
        T delta;
        std::vector<T> alpha;
        std::vector<T> beta;
        std::vector<T> C;
        std::vector<T> x;
        std::vector<T> grid;
        std::vector<T> grid_resp;
        std::vector<T> grid_wate;
        std::vector<T> P;
        std::vector<T> E;

    };

    ParksMcClellanParam *pmc_param;

    void setParksMcClellanParameters();
    void RemezInterchange(ParksMcClellanParam*);
    void ComputeLagrange(ParksMcClellanParam*);
    T EstimateFreqResponse(ParksMcClellanParam*, T);
    void ComputeCoeffs(ParksMcClellanParam*);

    T CalcFilterLen(const WindowType &w_type);
    void CalcIdealFIRCoeffs();
    void RectangularWindowCoeffs();
    void HammingWindowCoeffs();
    void HanningWindowCoeffs();
    void BlackmanWindowCoeffs();
    void KaiserWindowCoeffs(T beta);
    void BartlettWindowCoeffs();
    void ParksMcClellanWindowCoeffs();
    void MultWinIdealCoeffs();
public:

    DigFIRFilter(FilterType& f_type, FiltParam<T>& f_param, ApproxType& atype)
    {
        pmc_param = new ParksMcClellanParam();
    }

    ~DigFIRFilter()
    {
        delete pmc_param;
    }

    void CalcDigFIRCoeffs(const WindowType& w_type, std::size_t order = 0);

};

/*!
 * \brief
 * \param
 * \result
 */
template<typename T>
T DigFIRFilter<T>::CalcFilterLen(const WindowType& w_type)
{
    T result = 0.;
    T freq_delta, freq_delta_lower, freq_delta_upper;
    if(m_ftype == FilterType::LPF || m_ftype == FilterType::HPF)
    {
        freq_delta = std::fabs(m_fparam.freq_stopband.first()
                               - m_fparam.freq_passband.first()) / m_fparam.fsamp;
    }
    else if(m_ftype == FilterType::PBF || m_ftype == FilterType::PBF)
    {
        freq_delta_lower = std::fabs(m_fparam.freq_stopband.first()
                                     - m_fparam.freq_passband.first()) / m_fparam.fsamp;
        freq_delta_upper = std::fabs(m_fparam.freq_stopband.second()
                                     - m_fparam.freq_passband.second()) / m_fparam.fsamp;
        if(freq_delta_lower > freq_delta_upper)
        {
            freq_delta = freq_delta_upper;
        }
        else{
            freq_delta = freq_delta_lower;
        }
    }

    auto sb_error = std::pow(10.0, 0,05 * m_fparam.gain_stopband.first());
    auto pb_error = 1 - std::pow(10.0, 0,05 * m_fparam.gain_passband.first());

    T min_error;
    if(w_type == WindowType::RectangularWindow ||
       w_type == WindowType::HammingWindow ||
       w_type == WindowType::HanningWindow ||
       w_type == WindowType::BlackmanWindow ||
       w_type == WindowType::BartlettWindow ||
       w_type == WindowType::KaiserWindow)
    {
        if(sb_error < pb_error){
            min_error = sb_error;
        } else{
            min_error = pb_error;
        }

        auto db_error = -20 * std::log10(min_error);
        if(db_error > 50){
            result = 0.1102 * (db_error - 8.7);
        }
        else if(db_error > 21){
            result = 0.5842 * std::pow((db_error - 21), 0.4) + 0.07886 * (db_error - 21);
        }
        else{
            result = 0.;
        }

        if(db_error > 21){
            m_order = static_cast<std::size_t>(std::ceil((db_error - 7.95) / (2.285 * freq_delta))
                                               + 1);
        }
        else{
            m_order = static_cast<std::size_t>(std::ceil((5.794 / freq_delta) + 1));
        }

        if(!(m_order % 2)){
            ++m_order;
        }

    } else if(w_type == WindowType::ParksMcClellanWindow){
        sb_error = std::log10(sb_error);
        pb_error = std::log10(pb_error);
        freq_delta /= ADF_PI_2;

        auto k1 = (0.005309 * pb_error * pb_error + 0.07114 * pb_error - 0.4761) * sb_error
                - (0.00266 * pb_error * pb_error + 0.5941 * pb_error + 0.4278);
        auto k2 = 0.51244 * (pb_error - sb_error) = 11.012;
        m_order = static_cast<std::size_t>(std::ceil(((k1 - k2 * freq_delta * freq_delta)
                                                      / freq_delta) + 1.5));

        if(!(m_order % 2)){
            ++m_order;
        }
        result = 0.;
    }
    return result;
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::CalcIdealFIRCoeffs()
{
    auto tau = static_cast<T>(m_order  - 1) / 2;

    auto cut_freq1 = (m_fparam.freq_stopband.first() + m_fparam.freq_passband.first()) / (2 * m_fparam.fsamp);
    auto cut_freq2 = (m_fparam.freq_stopband.second() + m_fparam.freq_passband.second()) / (2 * m_fparam.fsamp);

    switch(m_ftype)
    {
    case FilterType::LPF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = cut_freq1 / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = std::sin(cut_freq1 * t) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::HPF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (ADF_PI - cut_freq1) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(ADF_PI * t) - std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::PBF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (cut_freq2 - cut_freq1) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(cut_freq2 * t) - std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;

    case FilterType::SBF :
        for(auto i=0; i<m_order; ++i)
        {
            if(i == static_cast<decltype(i)>(tau)){
                n_bcoefs[i] = (ADF_PI + cut_freq1 - cut_freq2) / ADF_PI;
            }
            else{
                auto t = i - tau;
                n_bcoefs[i] = (std::sin(ADF_PI * t) - std::sin(cut_freq2 * t)
                               + std::sin(cut_freq1 * t)) / (ADF_PI * t);
            }
        }
        break;
    }
}

template<typename T>
void DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType& w_type, std::size_t order)
{
    auto beta = CalcFilterLen(w_type);

    if(order < 0) return;
    if(order > 0)
    {
        m_order = order;
    }

    n_acoefs.reserve(m_order);
    n_bcoefs.reserve(m_order);
    n_acoefs.assign(m_order, 0);
    n_bcoefs.assign(m_order, 0);

    m_fparam.gain = 1.0;

    if(w_type != WindowType::ParksMcClellanWindow)
    {
        CalcIdealFIRCoeffs();
    }

    switch(w_type)
    {
    case WindowType::RectangularWindow :
        RectangularWindowCoeffs();
        break;

    case WindowType::BartlettWindow :
        BartlettWindowCoeffs();
        break;

    case WindowType::BlackmanWindow :
        BlackmanWindowCoeffs();
        break;

    case WindowType::HammingWindow :
        HammingWindowCoeffs();
        break;

    case WindowType::HanningWindow :
        HanningWindowCoeffs();
        break;

    case WindowType::KaiserWindow :
        KaiserWindowCoeffs(beta);
        break;

    case WindowType::ParksMcClellanWindow :
        ParksMcClellanWindowCoeffs();
        break;

    case WindowType::Undefined :
        ADF_ERROR("Error in DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType&, std::size_t): The unknown window type");
        return;

    default:
        ADF_ERROR("Error in DigFIRFilter<T>::CalcDigFIRCoeffs(const WindowType&, std::size_t): The window type value, is empty");
        return;
    }

    if(w_type != WindowType::ParksMcClellanWindow)
    {
        MultWinIdealCoeffs();
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::RectangularWindowCoeffs()
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = 1.0;
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::BartlettWindowCoeffs()
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = (2.0 * i) / (n - 1);
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::BlackmanWindowCoeffs()
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = 0.42 - 0.50 * std::cos(ADF_PI * i / (n - 1))
                              + 0.50 * std::cos(2 * ADF_PI * i / (n - 1));
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::HammingWindowCoeffs()
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = 0.54 - 0.46 * std::cos(ADF_PI * i / (n - 1));
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::HanningWindowCoeffs()
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = 0.50 - 0.50 * std::cos(ADF_PI * i / (n - 1));
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 * \param
 */
template<typename T>
void DigFIRFilter<T>::KaiserWindowCoeffs(T beta)
{
    auto N = m_order;
    auto n = static_cast<T>(N);
    for(auto indx=0; indx<(N+1)/2; ++indx)
    {
        auto i = static_cast<T>(indx);
        n_acoefs[indx] = bessel_func_mod(2 * beta * std::sqrt(i * (n - i -1)) / (n - 1)) /
                         bessel_func_mod(beta);
        n_acoefs[N - 1 - indx] = n_acoefs[indx];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::ParksMcClellanWindowCoeffs()
{
    pmc_param->grid_mult = ADF_GRID_MULT * 2;

    T min_err, max_err;
    for(auto k=0; k<4; ++k)
    {
        pmc_param->grid_mult *= 2;
        setParksMcClellanParameters();

        //! Start the Remez algorithm inner loop
        auto q_old = 1.0;
        bool done = false;
        for(auto i=0; i<30; ++i)
        {
            RemezInterchange(pmc_param);
            max_err = pmc_param->E[pmc_param->extrml[0]];
            min_err = max_err;
            for(auto j=1; j<pmc_param->num_ext; ++j)
            {
                auto test = pmc_param->E[pmc_param->extrml[j]];
                if(test > max_err)
                {
                    max_err = test;
                }
                if(test < min_err)
                {
                    min_err = test;
                }
            }
            auto q_now = static_cast<T>((max_err - min_err)) / max_err;
            if(q_now < ADF_ERR_EPS)
            {
                done = true;
                break;
            }

            //! If no significant progress after 15 iterations, try again with increase of grid density
            if((i > 15) && (q_now > q_old)){ break; }
            q_old = q_now;
        }

        if(done)
        {
            ComputeCoeffs(pmc_param);
            break;
        }
    }
}

template<typename T>
void DigFIRFilter<T>::MultWinIdealCoeffs()
{
    for(auto i=0; i<m_order; ++i)
    {
        n_acoefs[i] *= n_bcoefs[i];
    }
}

template<typename T>
void DigFIRFilter<T>::setParksMcClellanParameters()
{
    std::vector<T> edge((2*ADF_MAX_BANDS), 0);
    std::vector<T> resp(ADF_MAX_BANDS, 0);
    std::vector<T> wate(ADF_MAX_BANDS, 0);

    auto sb_err = std::pow(10, 0.05 * m_fparam.gain_stopband.first());
    auto pb_err = 1 - std::pow(10, 0.05 * m_fparam.gain_passband.first());
    auto W = pb_err / sb_err;

    switch(m_ftype)
    {
    case FilterType::LPF :
        pmc_param->num_bands = 2;
        edge[0] = 0;
        edge[1] = m_fparam.freq_passband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[2] = m_fparam.freq_stopband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[3] = 0.5;
        resp[0] = 1.0;
        resp[1] = 0.0;
        wate[0] = 1.0;
        wate[1] = W;
        break;

    case FilterType::HPF :
        pmc_param->num_bands = 2;
        edge[0] = 0;
        edge[1] = m_fparam.freq_stopband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[2] = m_fparam.freq_passband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[3] = 0.5;
        resp[0] = 0.0;
        resp[1] = 1.0;
        wate[0] = W;
        wate[1] = 1.0;
        break;

    case FilterType::PBF :
        pmc_param->num_bands = 3;
        edge[0] = 0;
        edge[1] = m_fparam.freq_stopband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[2] = m_fparam.freq_passband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[3] = m_fparam.freq_passband.second() / (ADF_PI_2 * m_fparam.fsamp);
        edge[4] = m_fparam.freq_stopband.second() / (ADF_PI_2 * m_fparam.fsamp);
        edge[5] = 0.5;
        resp[0] = 0.0;
        resp[1] = 1.0;
        resp[3] = 0.0;
        wate[0] = W;
        wate[1] = 1.0;
        wate[2] = W;
        break;

    case FilterType::SBF :
        pmc_param->num_bands = 3;
        edge[0] = 0;
        edge[1] = m_fparam.freq_passband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[2] = m_fparam.freq_stopband.first() / (ADF_PI_2 * m_fparam.fsamp);
        edge[3] = m_fparam.freq_stopband.second() / (ADF_PI_2 * m_fparam.fsamp);
        edge[4] = m_fparam.freq_passband.second() / (ADF_PI_2 * m_fparam.fsamp);
        edge[5] = 0.5;
        resp[0] = 1.0;
        resp[1] = 0.0;
        resp[3] = 1.0;
        wate[0] = 1.0;
        wate[1] = W;
        wate[2] = 1.0;
        break;
    }

    auto filt_len = m_order;
    auto num_cfs = (filt_len + 1) / 2;
    pmc_param->filt_type = FIRType::SymOdd;
    pmc_param->num_ext = num_cfs + 1;
    pmc_param->num_pts = pmc_param->num_ext * pmc_param->grid_mult;

    pmc_param->band_pts.reserve(ADF_MAX_BANDS);
    pmc_param->band_pts.assign(pmc_param->band_pts.capacity(), 0);

    pmc_param->grid.reserve(pmc_param->num_pts);
    pmc_param->grid.assign(pmc_param->grid.capacity(), 0);

    pmc_param->grid_resp.reserve(pmc_param->num_pts);
    pmc_param->grid_resp.assign(pmc_param->grid_resp.capacity(), 0);

    pmc_param->grid_wate.reserve(pmc_param->num_pts);
    pmc_param->grid_wate.assign(pmc_param->grid_wate.capacity(), 0);

    pmc_param->x.reserve(pmc_param->num_ext);
    pmc_param->x.assign(pmc_param->x.capacity(), 0);

    pmc_param->alpha.reserve(pmc_param->num_ext);
    pmc_param->alpha.assign(pmc_param->alpha.capacity(), 0);

    pmc_param->beta.reserve(pmc_param->num_ext);
    pmc_param->beta.assign(pmc_param->beta.capacity(), 0);

    pmc_param->C.reserve(pmc_param->num_ext);
    pmc_param->C.assign(pmc_param->C.capacity(), 0);

    pmc_param->extrml.reserve(2 * pmc_param->num_ext);
    pmc_param->extrml.assign(pmc_param->extrml.capacity(), 0);

    pmc_param->P.reserve(pmc_param->num_pts);
    pmc_param->P.assign(pmc_param->P.capacity(), 0);

    pmc_param->E.reserve(pmc_param->num_pts);
    pmc_param->E.assign(pmc_param->E.capacity(), 0);

    auto BW_total = 0;
    std::vector<T> band_wid(pmc_param->num_bands, 0);
    for(auto b=0; b<pmc_param->num_bands; ++b)
    {
        band_wid[b] = (edge[2*b + 1] - edge[2*b]);
        BW_total += band_wid[b];
    }

    auto grid_dens = static_cast<T>(pmc_param->num_pts - pmc_param->num_bands) / BW_total;
    auto extr_dens = static_cast<T>(pmc_param->num_ext - pmc_param->num_bands) / BW_total;
    std::vector<T> band_ext(ADF_MAX_BANDS, 0);
    band_ext[pmc_param->num_bands-1] = pmc_param->num_ext;
    pmc_param->band_pts[pmc_param->num_bands-1] = pmc_param->num_pts;

    for(auto b=0; b<pmc_param->num_bands-1; ++b)
    {
        pmc_param->band_pts[b] = static_cast<int>(grid_dens * band_wid[b] + 1.5);
        pmc_param->band_pts[pmc_param->num_bands-1] -= pmc_param->band_pts[b];
        band_ext[b] = static_cast<int>(extr_dens * band_wid[b] + 1.5);
        band_ext[pmc_param->num_bands-1] -= band_ext[b];
    }

    auto j=0;
    for(auto b=0; b<pmc_param->num_bands; ++b)
    {
        auto freq_spc = band_wid[b] / static_cast<T>(pmc_param->band_pts[b] - 1);
        pmc_param->grid[j] = edge[2 * b];
        pmc_param->grid_resp[j] = resp[b];
        pmc_param->grid_wate[j++] = wate[b];

        for(auto i=1; i<pmc_param->band_pts[b]; ++i)
        {
            pmc_param->grid[j] = pmc_param->grid[j-1] + freq_spc;
            pmc_param->grid_resp[j] = resp[b];
            pmc_param->grid_wate[j++] = wate[b];
        }
    }

    auto e = 0;
    auto index = 0;
    for(auto b=0; b<pmc_param->num_bands; ++b)
    {
        auto pt_dens = static_cast<T>(pmc_param->band_pts[b]) / static_cast<T>(band_ext[b]);
        pmc_param->extrml[e++] = index;

        for(auto i=0; i<band_ext[b]; ++i)
        {
            pmc_param->extrml[e++] = static_cast<int>(static_cast<T>(i) * pt_dens + 0.5) + index;
        }
        index += pmc_param->band_pts[b];
    }
}

/*!
 * \brief
 */
template<typename T>
void DigFIRFilter<T>::RemezInterchange(DigFIRFilter<T>::ParksMcClellanParam *param)
{
    //! Compute alpha[], beta[], C[], x[] and delta
    ComputeLagrange(param);

    for(auto j=0; j<param->num_pts; ++j)
    {
        param->P[j] = EstimateFreqResponse(param, ADF_PI_2 * param->grid[j]);
        param->E[j] = std::fabs((param->grid_resp[j] - param->P[j]) * param->grid_wate[j]);
    }

    T p, err_old, err_now;
    auto ext = 0;
    auto bnd_str = 0;
    for(auto b=0; b<param->num_pts; ++b)
    {
        p = EstimateFreqResponse(param, ADF_PI_2 * param->grid[bnd_str]);
        err_old = std::fabs((p - param->grid_resp[bnd_str]) * param->grid_wate[bnd_str]);

        p = EstimateFreqResponse(param, ADF_PI_2 * param->grid[bnd_str+1]);
        err_now = std::fabs((p - param->grid_resp[bnd_str+1]) * param->grid_wate[bnd_str+1]);

        auto slope = 1;
        if((b > 0) || (err_now < err_old))
        {
            param->extrml[ext++] = bnd_str;
            slope = -1;
        }

        err_old = err_now;
        auto j = 0;
        for(auto i=2; i<param->band_pts[b]-1; ++i)
        {
            j = i + bnd_str;
            p = EstimateFreqResponse(param, ADF_PI_2 * param->grid[j]);
            err_now = std::fabs((p - param->grid_resp[j]) * param->grid_wate[j]);

            if(err_now > err_old)
            {
                slope = -1;
            }
            else if(err_now < err_old)
            {
                if(slope == 1)
                {
                    param->extrml[ext++] = j - 1;
                    slope = -1;
                }
            }
            else
            {   }

            err_old = err_now;
        }

        p = EstimateFreqResponse(param, ADF_PI_2 * param->grid[j+1]);
        err_now = std::fabs((p - param->grid_resp[j+1]) * param->grid_wate[j+1]);

        if((b < param->num_bands-1) || (err_now > err_old))
        {
            param->extrml[ext++] = j+1;
        }
        bnd_str += param->band_pts[b];
    }

    auto excess = ext - param->num_ext;
    if(excess > 0)
    {
        for(auto i=0; i< excess; ++i)
        {
            auto index = 0;
            err_now = param->E[param->extrml[0]];

            for(auto j=1; j<ext; ++j)
            {
                if(param->E[param->extrml[j]] < err_now)
                {
                    index = j;
                    err_now = param->E[param->extrml[j]];
                }
            }
            --ext;

            for(auto j=index; j<ext; ++j)
            {
                param->extrml[j] = param->extrml[j+1];
            }
        }
    }
}

/*!
 * \brief ComputeLagrange - use Lagrange interpolation formula
 *                          \f$( A_{e}(e^{j\omega}) = P(\cos \omega) = \frac{\sum_{k=1}^{L+1}[d_{k}/(x - x_{k})]C_{k}}{\sum_{k=1}^{L+1}[d_{k}/(x - x_{k})]}
 *                          )\f$
 * \param param - struct of common values
 */
template<typename T>
void DigFIRFilter<T>::ComputeLagrange(DigFIRFilter<T>::ParksMcClellanParam *param)
{
    //! x_{i} = \cos(\omega_{i}), \omega_{i} are points of alternations
    for(auto k=0; k<param->num_ext; ++k)
    {
        param->x[k] = std::cos(param->grid[param->extrml[k]] * ADF_PI_2);
    }

    //! compute alpha[k] and beta[k] see b_{k}, d_{k} values in Oppenheim A.V.-Discrete-Time Signal Processing.-2013 pp.591
    for(auto k=0; k<param->num_ext; ++k)
    {
        auto prod = 1.0;
        auto xk = param->x[k];
        for(auto i=0; i<param->num_ext; ++i)
        {
            if(xk != param->x[i]){
                prod *= (xk - param->x[i]);
            }
        }
        if(prod == 0.0){
            param->alpha[k] = ADF_INF;
        }
        else{
            param->alpha[k] = 1.0/prod;
        }

        param->beta[k] = param->alpha[k] * (xk - param->x[param->num_ext-1]);
    }

    //! compute \delta see Oppenheim A.V.-Discrete-Time Signal Processing.-2013 pp.591 eq.114
    int index;
    auto numer = 0.0;
    auto denom = 0.0;
    auto sign = 1.0;
    for(auto k=0; k<param->num_ext; ++k)
    {
        index = param->extrml[k];
        numer += param->alpha[k] * param->grid_resp[index];
        denom += sign * param->alpha[k] / param->grid_wate[index];
        sign = -sign;
    }

    if(denom == 0.0){
        param->delta = ADF_INF;
    }
    else{
        param->delta = numer / denom;
    }
    sign = 1.0;

    //! compute C[k} in Oppenheim A.V.-Discrete-Time Signal Processing.-2013 pp.591 eq.116b
    for(auto k=0; k<param->num_ext; ++k)
    {
        index = param->extrml[k];
        param->C[k] = param->grid_resp[index] - sign * param->delta /
                                               param->grid_wate[index];
        sign = -sign;
    }
}

template<typename T>
T DigFIRFilter<T>::EstimateFreqResponse(DigFIRFilter::ParksMcClellanParam *param, T rad_freq)
{
    T denom_sum, numer_sum;
    auto xf = std::cos(rad_freq);
    for(auto k=0; k<param->num_ext; ++k)
    {
        if(std::fabs(xf - param->x[k] > ERR_SMALL))
        {
            auto tmp = param->beta[k] / (xf - param->x[k]);
            denom_sum += tmp;
            numer_sum += tmp * param->C[k];
        }
        else{
            return param->C[k];
        }
    }

    if(denom_sum == 0.0)
    {
        return static_cast<T>(ADF_INF);
    }
    else{
        return static_cast<T>(numer_sum/denom_sum);
    }
}

template<typename T>
void DigFIRFilter<T>::ComputeCoeffs(DigFIRFilter::ParksMcClellanParam *param)
{
    auto order = m_order;
    auto num_cfs = (order + 1) / 2;
    auto omega = ADF_PI_2 / static_cast<T>(order);
    for(auto k=0; k<num_cfs; ++k)
    {
        param->P[k] = EstimateFreqResponse(param, k*omega);
    }

    switch(param->filt_type)
    {
    case FIRType::SymOdd :
        for(auto n=0; n<num_cfs; ++n)
        {
            auto m = num_cfs-1-n;
            n_acoefs[m] = param->P[0];
            auto mult = ADF_PI_2 * n / static_cast<T>(order);

            for(auto k=1; k<num_cfs; ++k)
            {
                n_acoefs[m] += 2 * param->P[k] * std::cos(k * mult);
            }
            n_acoefs[m] /= static_cast<T>(order);
            n_acoefs[order-1-m] = n_acoefs[m];
        }
        break;
    case FIRType::SymEven :
        break;
    case FIRType::AsymOdd :
        break;
    case FIRType::AsymEven :
        break;
    default:
        ADF_ERROR("Error in DigFIRFilter<T>::ComputeCoeffs(DigFIRFilter::ParksMcClellanParam *) Unknown FIR filter type");
    }
}

} //namespace adf
#endif //DIGFIRFILTER_H
