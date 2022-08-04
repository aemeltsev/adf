// Copyright (C) 2021-2022 Anthony Emeltsev
// SPDX-License-Identifier: Apache-2.0
//

#ifndef FFT_H
#define FFT_H

#include <cmath>
#include <cstdint>
#include "base.hpp"
#include "complex.hpp"

namespace adf{
/**
 * \brief FFT radix-r 1 dimension non-recursive algorithm
 */
template<typename T = double>
class FFT
{
private:
    int m_vec_size;
    int m_log2_vec_size;
    enum m_direction {forward, inverse};
    complex<T> *m_uforward, *m_uinverse;
    complex<T> *m_workspace;
    int *m_permindex;

    /*!
     * \brief swap - Helper method for swapping two complex values
     * \param a - complex value
     * \param b - complex value
     */
    inline void swap(complex<T> &a, complex<T> &b)
    {
        std::swap(a, b);
    }

    /*!
     * \brief Setup - Helper method for checking that the size of the array
     *                of samples is not zero and equal to the power of 2.
     *                Also, the initialization of the array of rotation factors
     *                and the array of indices for permutations
     * \param size - size of samples array
     */
    void Setup(int size)
    {
        if(size == m_vec_size) return;
        if(size < 1) ADF_ERROR("Error in FFT::Setup(int): Requested length must be > 0");

        int k;
        //! Check that size is power of 2.
        for(k = size, m_log2_vec_size = 1; k > 2; k /= 2, ++m_log2_vec_size)
        {
            if((k % 2) != 0) ADF_ERROR("Error in FFT::Setup(int): Requested length is not a power of 2");
        }

        ReleaseMemory();
        m_vec_size = size;

        //! Allocate and setup complex arrays.
        if(!(m_uforward = new complex<T>[size])) ADF_ERROR("Error in FFT::Setup: ErrNoMem");
        if(!(m_uinverse = new complex<T>[size])) ADF_ERROR("Error in FFT::Setup: ErrNoMem");;
        if(!(m_workspace = new complex<T>[size])) ADF_ERROR("Error in FFT::Setup: ErrNoMem");;

        T baseang = -2 * ADF_PI / static_cast<T>(size);

        m_uforward[0] = m_uinverse[0] = complex<T>(1,0);

        for(k = 1; k < size / 2; ++k)
        {
            auto x = std::cos(baseang * k);
            auto y = -std::sqrt(1 - x * x);
            m_uforward[k] = m_uinverse[size - k] = complex<T>(x, y);
            m_uforward[size - k] = m_uinverse[k] = complex<T>(x, -y);
        }

        if(size > 1)
        {
            m_uforward[size / 2] = m_uinverse[size / 2] = complex<T>(-1, 0);
        }

        //! Allocate and setup(bit-reversal) permutation index.
        if(!(m_permindex = new int[size])) ADF_ERROR("Error in FFT::Setup ErrNoMem");

        m_permindex[0] = 0;
        int m, n;    /*!< The following code relies heavily on size == 2 ^ m_log2_vec_size */
        for(k = 1, n = (size >> 1); k < size; ++k)
        {
            //! At each step, n is bit-reversed pattern of k.
            if(n > k) m_permindex[k] = n;  /*!< Set a swap index */
            else m_permindex[k] = 0;       /*!< Do nothing: Index already swapped or the same */

            //! Calculate next n.
            m = (size >> 1);
            while((m > 0) && (n & m))
            {
                n -= m;
                m >>= 1;
            }
            n += m;
        }
    }

    /*!
     * \brief Permute - initial permutation
     * \param vec - initial array of the samples
     */
    void Permute(complex<T> *vec)
    {
        for(auto i = 0; i < m_vec_size; ++i)
        {
            int j;
            if((j = m_permindex[i]) != 0)
                swap(vec[i], vec[j]);
        }
    }

    /*!
     * \brief BaseDIFForward - In-place forward FFT using Decimation in Frequency technique,
     *                         without reshuffling of indices.
     *                         This code does not use the complex<T> class,
     *                         because operations on complex numbers such as multiplication
     *                         can be optimized not efficiently.
     *                         So this routine just makes use of ordinary type "double" variables,
     *                         and assumes each Complex variable is actually two consecutive double variables.

     * \param vec - initial array of the samples
     */
    void BaseDIFForward(complex<T> *vec)
    {
        if(m_vec_size == 1) return; /*!< Nothing to do */

        double *v;

        // TODO maybe need creating of a helper function for
        // convert raw array of the complex values to raw array T values
        double const *const U = static_cast<double *>(m_uforward);
        double *const dvec = static_cast<double *>(vec);

        int block, blocksize, blockcount, offset, uoff1;
        int halfbs, threehalfbs; /*!< Half blocksize, 3/2 blocksize */
        double m1x, m1y, m2x, m2y, m3x, m3y;
        double x0, y0, x1, y1, x2, y2, x3, y3;
        double xs02, ys02, xd02, yd02, xs13, ys13, xd13, yd13;
        double t1x, t1y, t2x, t2y, t3x, t3y;

        //! Blocksize > 4.
        for(blocksize = m_vec_size, blockcount = 1;
            blocksize > 4;
            blocksize /= 4, blockcount *= 4)
        {
            //! Loop through double-step matrix multiplications.
            halfbs = blocksize / 2;
            threehalfbs = blocksize + halfbs;

            for(block = 0, v = dvec; block < blockcount; ++block, v += (2 * blocksize))
            {
                for(offset = 0; offset < halfbs; offset += 2)
                {
                    uoff1 = offset * blockcount;
                    m1x = U[uoff1];
                    m1y = U[uoff1 + 1];

                    m2x = U[2 * uoff1];
                    m2y = U[2 * uoff1 + 1];

                    m3x = U[3 * uoff1];
                    m3y = U[3 * uoff1 + 1];

                    x0 = v[offset];
                    y0 = v[offset + 1];
                    x1 = v[halfbs + offset];
                    y1 = v[halfbs + offset + 1];
                    x2 = v[blocksize + offset];
                    y2 = v[blocksize + offset + 1];
                    x3 = v[threehalfbs + offset];
                    y3 = v[threehalfbs + offset + 1];

                    xs02 = x0 + x2;
                    xs13 = x1 + x3;
                    v[offset] = xs02 + xs13;
                    t1x = xs02 - xs13;

                    ys02 = y0 + y2;
                    ys13 = y1 + y3;
                    v[offset + 1] = ys02 + ys13;
                    t1y = ys02 - ys13;

                    v[halfbs + offset]      = t1x * m2x - t1y * m2y;
                    v[halfbs + offset + 1]  = t1y * m2x + t1x * m2y;

                    xd02 = x0 - x2;
                    yd13 = y1 - y3;
                    t3x = xd02 - yd13;
                    t2x = xd02 + yd13;
                    yd02 = y0 - y2;
                    xd13 = x1 - x3;
                    t3y = yd02 + xd13;
                    t2y = yd02 - xd13;
                    v[blocksize + offset]         = t2x * m1x - t2y * m1y;
                    v[blocksize + offset + 1]     = t2y * m1x + t2x * m1y;
                    v[threehalfbs + offset]       = t3x * m3x - t3y * m3y;
                    v[threehalfbs + offset + 1]   = t3y * m3x + t3x * m3y;
                }
            }
        }

        //! Do smallest blocks; size is either 4 or 2.
        if(blocksize == 4)
        {
            blockcount = m_vec_size / 4;
            for(block = 0, v = dvec; block < blockcount; ++block, v += 8)
            {
                x0 = v[0];
                y0 = v[1];

                x1 = v[2];
                y1 = v[3];

                x2 = v[4];
                y2 = v[5];

                x3 = v[6];
                y3 = v[7];

                xs02 = x0 + x2;
                xs13 = x1 + x3;
                v[0] = xs02 + xs13;
                v[2] = xs02 - xs13;

                ys02 = y0 + y2;
                ys13 = y1 + y3;
                v[1] = ys02 + ys13;
                v[3] = ys02 - ys13;

                xd02 = x0 - x2;
                yd13 = y1 - y3;
                v[4] = xd02 + yd13;
                v[6] = xd02 - yd13;

                yd02 = y0 - y2;
                xd13 = x1 - x3;
                v[5] = yd02 - xd13;
                v[7] = yd02 + xd13;
            }
        }
        else
        { //! blocksize == 2.
            blockcount = m_vec_size / 2;
            for(block = 0, v = dvec; block < blockcount; ++block, v += 4)
            {
                x0 = v[0];
                y0 = v[1];

                x1 = v[2];
                y1 = v[3];

                v[0] = x0 + x1;
                v[2] = x0 - x1;

                v[1] = y0 + y1;
                v[3] = y0 - y1;
            }
        }
    }

    /*!
     * \brief BaseDITInverse - In-place inverse FFT using Decimation in Time technique,
     *                         without reshuffling of indices.
     * \param vec - initial array of the samples
     */
    void BaseDITInverse(complex<T> *vec)
    {
        if(m_vec_size == 1) return; /*!< Nothing to do */

        double *v;

        // TODO maybe need creating of a helper function for
        // convert raw array of the complex values to raw array T values
        double const *const U = static_cast<double *>(m_uforward);
        double *const dvec = static_cast<double *>(vec);

        int block, blocksize, blockcount, offset, uoff1;
        int halfbs, threehalfbs; /*!< Half blocksize,3/2 blocksize */
        double m1x, m1y, m2x, m2y, m3x, m3y;
        double x0, y0, x1, y1, x2, y2, x3, y3;
        double xs01, ys01, xd01, yd01, xs23, ys23, xd23, yd23;
        double t1x, t1y, t2x, t2y, t3x, t3y;

        //! Do smallest blocks; size is either 4 or 2.
        if(m_vec_size > 2)
        {
            blockcount = m_vec_size / 4;
            for(block = 0, v = dvec; block < blockcount; ++block, v += 8)
            {
                x0 = v[0];
                y0 = v[1];

                x1 = v[2];
                y1 = v[3];

                x2 = v[4];
                y2 = v[5];

                x3 = v[6];
                y3 = v[7];

                xs01 = x0 + x1;
                xs23 = x2 + x3;
                v[0] = xs01 + xs23;
                v[4] = xs01 - xs23;

                ys01 = y0 + y1;
                ys23 = y2 + y3;
                v[1] = ys01 + ys23;
                v[5] = ys01 - ys23;

                xd01 = x0 - x1;
                yd23 = y2 - y3;
                v[2] = xd01 - yd23;
                v[6] = xd01 + yd23;

                yd01 = y0 - y1;
                xd23 = x2 - x3;
                v[3] = yd01 + xd23;
                v[7] = yd01 - xd23;
            }
        }
        else
        {//! vecsize == 2.
            x0 = dvec[0];
            y0 = dvec[1];

            x1 = dvec[2];
            y1 = dvec[3];

            dvec[0] = x0 + x1;
            dvec[2] = x0 - x1;
            dvec[1] = y0 + y1;
            dvec[3] = y0 - y1;
            return;
        }

        //! Blocksize > 4.
        for(blocksize = 16, blockcount = m_vec_size / 16;
            blocksize <= m_vec_size;
            blocksize *= 4, blockcount /= 4)
        {
            //! Loop through double-step matric multiplications.
            halfbs = blocksize / 2;
            threehalfbs = blocksize + halfbs;
            for(block = 0, v = dvec; block < blockcount; ++block, v += (2 * blocksize))
            {
                for(offset = 0; offset < blocksize / 2; offset += 2)
                {
                    x0 = v[offset];
                    y0 = v[offset + 1];
                    t2x = v[blocksize + offset];
                    t2y = v[blocksize + offset + 1];
                    uoff1 = offset * blockcount;
                    m1x = U[uoff1];
                    m1y = U[uoff1 + 1];
                    x2 = t2x * m1x - t2y * m1y;
                    y2 = t2y * m1x + t2x * m1y;

                    m2x = U[2 * uoff1];
                    m2y = U[2 * uoff1 + 1];
                    t1x = v[halfbs + offset];
                    t1y = v[halfbs + offset + 1];
                    x1 = t1x * m2x - t1y * m2y;
                    y1 = t1y * m2x + t1x * m2y;

                    t3x = v[threehalfbs + offset];
                    t3y = v[threehalfbs + offset + 1];
                    m3x = U[3 * uoff1];
                    m3y = U[3 * uoff1 + 1];
                    x3 = t3x * m3x - t3y * m3y;
                    y3 = t3y * m3x + t3x * m3y;

                    xs01 = x0 + x1;
                    xs23 = x2 + x3;
                    v[offset] = xs01 + xs23;
                    v[blocksize + offset] = xs01 - xs23;

                    ys01 = y0 + y1;
                    ys23 = y2 + y3;
                    v[offset + 1] = ys01 + ys23;
                    v[blocksize + offset + 1] = ys01 - ys23;

                    yd01 = y0 - y1;
                    xd23 = x2 - x3;
                    v[halfbs + offset + 1] = yd01 + xd23;
                    v[threehalfbs + offset + 1] = yd01 - xd23;

                    xd01 = x0 - x1;
                    yd23 = y2 - y3;
                    v[halfbs + offset] = xd01 - yd23;
                    v[threehalfbs + offset] = xd01 + yd23;
                }
            }
        }

        if(blocksize == 2 * m_vec_size)
        {
            //! We still have to do one single-step matrix multiplication.
            blocksize = m_vec_size;
            v = dvec;
            for(offset = 0; offset < m_vec_size; offset += 2, v += 2)
            {
                x0 = v[0];
                y0 = v[1];
                t1x = v[m_vec_size];
                t1y = v[m_vec_size + 1];
                m1x = U[offset];
                m1y = U[offset + 1];
                x1 = t1x * m1x - t1y * m1y;
                y1 = t1y * m1x + t1x * m1y;
                v[0] = x0 + x1;
                v[m_vec_size] = x0 - x1;
                v[1] = y0 + y1;
                v[m_vec_size + 1] = y0 - y1;
            }
        }
    }

public:

    /*!
     * \brief ctor
     */
    FFT()
    {
        m_vec_size = 0;
        m_uforward = m_uinverse = m_workspace = nullptr;
        m_permindex = nullptr;
    }

    /*!
     * \brief dtor using helper method for release memory
     */
    ~FFT()
    {
        ReleaseMemory();
    }

    /*!
     * \brief ForwardDecFreq - Forward decimation in frequency fft algorithm
     * \param size - size of samples array
     * \param vec - initial array of the samples
     * \param divisor - so to normalise (in any sense) the output data you need to multiply it by 1/N
     *
     */
    void ForwardDecFreq(int size, complex<T> *vec, T divisor=0.)
    {
        if(divisor == 0) divisor = 1.; /*!< Default is no normalization on forward FFT. */
        Setup(size);
        BaseDIFForward(vec);
        Permute(vec);
        if(divisor != 0. && divisor != 1.)
        {
            T mult = 1. / divisor;
            for(int k = 0; k < size; ++k) vec[k] *= mult;
        }
    }

    /*!
     * \brief InverseDecTime - Inverse decimation in time fft algorithm
     * \param size - size of samples array
     * \param vec - initial array of the samples
     * \param divisor - so to normalise (in any sense) the output data you need to multiply it by 1/N
     */
    void InverseDecTime(int size, complex<T> *vec, T divisor=0.)
    {
        if(divisor == 0) divisor = static_cast<T>(size); /*!< Default divisor on iFFT is 'size'. */
        Setup(size);
        Permute(vec);
        BaseDITInverse(vec);
        if(divisor != 0. && divisor != 1.)
        {
            T mult = 1. / divisor;
            for(int k = 0; k < size; ++k) vec[k] *= mult;
        }
    }

    /*!
     * \brief ReleaseMemory - clean allocated memory
     */
    void ReleaseMemory()
    {
        if(m_vec_size > 0)
        {
            if(m_uforward != nullptr) delete[] m_uforward;
            if(m_uinverse != nullptr) delete[] m_uinverse;
            if(m_workspace != nullptr) delete[] m_workspace;
            if(m_permindex != nullptr) delete[] m_permindex;
        }

        m_vec_size = 0;
        m_uforward = m_uinverse = m_workspace = nullptr;
        m_permindex = nullptr;
    }
};

/*!
 *
 */
template<typename T = double>
class FFTR2D
{
private:
    // Relative speed of ForwardCR as compared to ForwardRC.
    // If bigger than 1, then ForwardCR is
    // faster.  This will be machine & compiler dependent...oh well.
    static constexpr double m_crrc_speed_ratio = 1.10;

    int m_vec_size1, m_vec_size2; /*!< Dimensions of full-sized real array.*/
    int m_log_size1, m_log_size2; /*!< Base-2 logs of m_vec_size1 & m_vec_size2 */
                                  /*!< (which *must* be powers of 2).*/
                                  /*!< The corresponding MyComplex array is only half-height, i.e., */
                                  /*!< is (m_vec_size1 / 2) + 1 rows by m_vec_size2 columns. */
    FFT<T> m_fft1, m_fft2;
    complex<T> *m_scratcha;
    complex<T> *m_scratchb;  /*!< Secondary scratch space; same size as 'scratch'. */
                             /*!< The main use of m_scratchb is to slightly block memory accesses */
                             /*!< having bad strides. */
    complex<T> **m_work_arr; /*!< Used only by inverse FFT routines */

public:

    /*!
     * \brief ReleaseMemory
     */
    void ReleaseMemory()
    {
        if(m_vec_size1 == 0 || m_vec_size2 == 0) return;

        delete[] m_scratcha;
        m_scratcha = nullptr;

        delete[] m_scratchb;
        m_scratchb = nullptr;

        if(m_work_arr != nullptr)
        {
            delete[] m_work_arr[0];
            delete[] m_work_arr;
            m_work_arr = nullptr;
        }

        m_vec_size1 = 0;
        m_vec_size2 = 0;
        m_fft1.ReleaseMemory(); /*!< m_fft1 and m_fft2 stack variable which inside uses a raw pointers*/
        m_fft2.ReleaseMemory();
    }

private:

    /*!
     * \brief Setup
     * \param size1
     * \param size2
     */
    void Setup(int size1, int size2)
    {
        if(size1 == m_vec_size1 && size2 == m_vec_size2) return;  /*!< Nothing to do */

        //! Release memory on size = 0 request.
        if(size1 == 0 || size2 == 0)
        {
            ReleaseMemory();
            return;
        }

        //! Check that sizes are powers of 2, and > 1. Also extract base-2 log of sizes.
        if(size1 < 2 || size2 < 2) ADF_ERROR("Error in FFTReal2D::Setup(int): Requested size1 or size2 must be >=2");

        int k;
        for(k = size1, m_log_size1 = 1; k > 2; k /= 2, ++m_log_size1)
        {
            if((k % 2) != 0)
                ADF_ERROR("Error in FFTReal2D::Setup(int): Requested size1 is not a power of 2");
        }
        for(k = size2, m_log_size2 = 1; k > 2; k /= 2, ++m_log_size2)
        {
            if((k % 2) != 0)
                ADF_ERROR("Error in FFTReal2D::Setup(int): Requested size2 is not a power of 2");
        }

        //! Allocate new space.
        ReleaseMemory();
        m_scratcha = new complex<T>[std::max(size1, size2)];
        m_scratchb = new complex<T>[std::max(size1, size2)];
        m_vec_size1 = size1;
        m_vec_size2 = size2;
    }

    /*!
     * \brief SetupInverse
     * \param size1
     * \param size2
     */
    void SetupInverse(int size1, int size2)
    {
        if(size1 == m_vec_size1 && size2 == m_vec_size2 && m_work_arr != nullptr)
            return;  /*!< Nothing to do */

        Setup(size1,size2);

        if(m_work_arr != nullptr)
        {
            delete[] m_work_arr[0];
            delete[] m_work_arr;
        } // Safety

        auto rowcount = (m_vec_size1 / 2) + 1;
        m_work_arr = new complex<T> * [rowcount];
        m_work_arr[0] = new complex<T>[rowcount * m_vec_size2];
        for(auto i = 1; i < rowcount; ++i)
        {
            m_work_arr[i] = m_work_arr[i - 1] + m_vec_size2;
        }
    }

    /*!
     * \brief ForwardCR - Routine ForwardCR does an FFT first on the columns,
     *                    then on the rows. It is natural to pair ForwardCR with InverseRC.
     *                    The routine choice is made automatically
     *                    (based on speed estimates)by the routines ::Forward &
     *                    ::Inverse (below in the "public" access block).
     * \param rsize1 - dimension of input array
     * \param rsize2 - dimension of input array
     * \param rarr - input 2 dimension array
     * \param csize1 - dimension of output array
     * \param csize2 - dimension of output array
     * \param carr - output result array with complex values
     */
    void ForwardCR(int rsize1,
                   int rsize2,
                   const double* const* rarr,
                   int csize1,
                   int csize2,
                   complex<T>** carr)
    {
        Setup(2 * (csize1 - 1), csize2); /*!< Safety */
        if(m_vec_size1 == 0 || m_vec_size2 == 0) return; /*!< Nothing to do */

        int i, j;
        double x1, x2, y1, y2;
        double xb1, xb2, yb1, yb2;

        //! Do FFT on columns, 2 at a time.
        for(j = 0; j + 3 < rsize2; j += 4)
        {
            //! Pack into complex<T> scratch array.
            for(i = 0; i < rsize1; ++i)
            {
                m_scratcha[i] = complex<T>(rarr[i][j],     rarr[i][j + 1]);
                m_scratchb[i] = complex<T>(rarr[i][j + 2], rarr[i][j + 3]);
            }

            //! Zero pad scratch space.
            for(i = rsize1; i < m_vec_size1; ++i) m_scratcha[i] = complex<T>(0., 0.);
            for(i = rsize1; i < m_vec_size1; ++i) m_scratchb[i] = complex<T>(0., 0.);

            //! Do complex FFT.
            m_fft1.ForwardDecFreq(m_vec_size1, m_scratchb);
            m_fft1.ForwardDecFreq(m_vec_size1, m_scratcha);

            //! Unpack into top half of 2D complex array.
            //! Rows 0 & m_vec_size1 / 2 are real-valued, so pack them together
            //! into row 0 (row 0 as real part, row m_vec_size1 / 2 as imag. part).
            carr[0][j    ] = complex<T>(m_scratcha[0].getReal(), m_scratcha[m_vec_size1 / 2].getReal());
            carr[0][j + 1] = complex<T>(m_scratcha[0].getImag(), m_scratcha[m_vec_size1 / 2].getImag());
            carr[0][j + 2] = complex<T>(m_scratchb[0].getReal(), m_scratchb[m_vec_size1 / 2].getReal());
            carr[0][j + 3] = complex<T>(m_scratchb[0].getImag(), m_scratchb[m_vec_size1 / 2].getImag());

            for(i = 1; i < m_vec_size1 / 2; ++i)
            { //! ASSUMES vecsize1 is even!
                x1 = m_scratcha[i].getReal() / 2;
                y1 = m_scratcha[i].getImag() / 2;

                x2 = m_scratcha[m_vec_size1 - i].getReal() / 2;
                y2 = m_scratcha[m_vec_size1 - i].getImag() / 2;

                xb1 = m_scratchb[i].getReal() / 2;
                yb1 = m_scratchb[i].getImag() / 2;

                xb2 = m_scratchb[m_vec_size1 - i].getReal() / 2;
                yb2 = m_scratchb[m_vec_size1 - i].getImag() / 2;

                carr[i][j    ] = complex<T>(x1 + x2,    y1 - y2);
                carr[i][j + 1] = complex<T>(y1 + y2,    x2 - x1);
                carr[i][j + 2] = complex<T>(xb1 + xb2, yb1 - yb2);
                carr[i][j + 3] = complex<T>(yb1 + yb2, xb2 - xb1);
            }
        }

        //! Case rsize2 not divisible by 4.
        for(; j < rsize2; j += 2)
        {
            //! Pack into complex scratch array.
            if(j + 1 < rsize2)
            {
                for(i = 0; i < rsize1; ++i)
                {
                    m_scratcha[i] = complex<T>(rarr[i][j], rarr[i][j + 1]);
                }
            }
            else { //! rsize2 == 1 mod 2.
                for(i = 0; i < rsize1; ++i)
                {
                    m_scratcha[i] = complex<T>(rarr[i][j], 0.);
                }
            }

            for(i = rsize1; i < m_vec_size1; ++i) m_scratcha[i] = complex<T>(0., 0.);

            m_fft1.ForwardDecFreq(m_vec_size1, m_scratcha);
            carr[0][j    ] = complex<T>(m_scratcha[0].getReal(),
                                        m_scratcha[m_vec_size1 / 2].getReal());

            carr[0][j + 1] = complex<T>(m_scratcha[0].getImag(),
                                        m_scratcha[m_vec_size1/2].getImag());

            for(i = 1; i < m_vec_size1 / 2; ++i)
            { //! ASSUMES vecsize1 is even!
                x1 = m_scratcha[i].getReal() / 2;
                y1 = m_scratcha[i].getImag() / 2;

                x2 = m_scratcha[m_vec_size1 - i].getReal() / 2;
                y2 = m_scratcha[m_vec_size1 - i].getImag() / 2;

                carr[i][j    ] = complex<T>(x1 + x2, y1 - y2);
                carr[i][j + 1] = complex<T>(y1 + y2, x2 - x1);
            }
        }

        //! Zero-pad remaining columns.
        if(rsize2 < csize2)
        {
            for(i = 0; i < csize1; ++i) {
                for(j = rsize2; j < csize2; ++j) {
                    carr[i][j] = complex<T>(0., 0.);
                }
            }
        }

        //! Do FFT on top half of rows.
        for(i = 0; i < csize1 - 1; ++i) m_fft2.ForwardDecFreq(csize2, carr[i]);

        //! Pull out row 0 & row csize1-1 from (packed) row 0.
        carr[csize1 - 1][0] = complex<T>(carr[0][0].getImag(), 0.);
        carr[0][0]          = complex<T>(carr[0][0].getReal(), 0.);

        for(j = 1; j < csize2 / 2; ++j)
        {
            x1 = carr[0][j].getReal() / 2;
            y1 = carr[0][j].getImag() / 2;

            x2 = carr[0][csize2 - j].getReal() / 2;
            y2 = carr[0][csize2 - j].getImag() / 2;

            complex<T> temp1(x1 + x2, y1 - y2);
            complex<T> temp2(y1 + y2, x2 - x1);

            carr[0][j]                    = temp1;
            carr[0][csize2 - j]           = temp1.conj();
            carr[csize1 - 1][j]           = temp2;
            carr[csize1 - 1][csize2 - j]  = temp2.conj();
        }

        carr[csize1 - 1][csize2 / 2] = complex<T>(carr[0][csize2 / 2].getImag(), 0.);
        carr[0][csize2 / 2]          = complex<T>(carr[0][csize2 / 2].getReal(), 0.);
    }

    /*!
     * \brief ForwardRC - Routine ForwardRC does an FFT first on the rows,
     *                    and then on the columns. It is natural to pair ForwardRC with InverseCR
     *                    The routine choice is made automatically
     *                    (based on speed estimates)by the routines ::Forward &
     *                    ::Inverse (below in the "public" access block).
     * \param rsize1 - dimension of input array
     * \param rsize2 - dimension of input array
     * \param rarr - input 2 dimension array
     * \param csize1 - dimension of output array
     * \param csize2 - dimension of output array
     * \param carr - output result array with complex values
     */
    void ForwardRC(int rsize1,
                   int rsize2,
                   const double* const* rarr,
                   int csize1,
                   int csize2,
                   complex<T>** carr)
    {
        Setup(2 * (csize1 - 1), csize2); /*!< Safety */
        if(m_vec_size1 == 0 || m_vec_size2 == 0) return; /*!< Nothing to do */

        int i, j;
        double x1, x2, y1, y2;
        double xb1, xb2, yb1, yb2;

        //! Do row FFT's.
        for(i = 0; i + 1 < rsize1; i += 2)
        {
            //! Pack complex<T> row.
          for(j = 0; j < rsize2; ++j)
          {
              carr[i / 2][j] = complex<T>(rarr[i][j], rarr[i + 1][j]);
          }

          for(j = rsize2; j < csize2; ++j)
          {
              carr[i / 2][j] = complex<T>(0.,0.); /*!< Zero pad */
          }
          //! Do FFT.
          m_fft2.ForwardDecFreq(csize2, carr[i / 2]);
        }

        for(; i < rsize1; i += 2)
        { //! In case rsize1 == 1 mod 2.
          for(j = 0; j < rsize2; ++j)
          {
              carr[i / 2][j] = complex<T>(rarr[i][j], 0.);
          }

          for(j = rsize2; j < csize2; ++j)
          {
              carr[i / 2][j] = complex<T>(0., 0.); /*!< Zero pad */
          }
          //! Do FFT.
          m_fft2.ForwardDecFreq(csize2, carr[i / 2]);
        }
        //! Any remaining rows are zero padding on the fly during the column FFT's (see below).

        //! Do column FFT's.
        //! Do column 0 and csize2/2, making use of the fact that these 2 columns are 'real'.
        for(i = 0; i < (rsize1 + 1) / 2; ++i)
        {
            x1 = carr[i][0].getReal();
            x2 = carr[i][0].getImag();

            y1 = carr[i][csize2 / 2].getReal();
            y2 = carr[i][csize2 / 2].getImag();

            m_scratcha[2 * i      ] = complex<T>(x1, y1);
            m_scratcha[(2 * i) + 1] = complex<T>(x2, y2);
        }

        for(i *= 2; i < m_vec_size1; ++i)
        {
            m_scratcha[i] = complex<T>(0., 0.); /*!< Zero pad */
        }

        m_fft1.ForwardDecFreq(m_vec_size1, m_scratcha);

        carr[0][0         ] = complex<T>(m_scratcha[0].getReal(), 0.);
        carr[0][csize2 / 2] = complex<T>(m_scratcha[0].getImag(), 0.);

        for(i = 1; i < csize1 - 1; ++i)
        {
          x1 = m_scratcha[i].getReal() / 2;
          y1 = m_scratcha[i].getImag() / 2;

          x2 = m_scratcha[m_vec_size1 - i].getReal() / 2;
          y2 = m_scratcha[m_vec_size1 - i].getImag() / 2;

          carr[i][0         ] = complex<T>(x1 + x2, y1 - y2);
          carr[i][csize2 / 2] = complex<T>(y1 + y2, x2 - x1);
        }

        carr[csize1 - 1][0         ] = complex<T>(m_scratcha[csize1 - 1].getReal(), 0.);
        carr[csize1 - 1][csize2 / 2] = complex<T>(m_scratcha[csize1 - 1].getImag(), 0.);

        //! Do remaining columns.
        for(j = 1; j + 1 < csize2 / 2; j += 2)
        {
            for(i = 0; i < (rsize1 + 1) / 2; ++i)
            {
                x1 = carr[i][j].getReal() / 2;
                y1 = carr[i][j].getImag() / 2;

                xb1 = carr[i][j + 1].getReal() / 2;
                yb1 = carr[i][j + 1].getImag() / 2;

                xb2 = carr[i][csize2 - 1 - j].getReal() / 2;
                yb2 = carr[i][csize2 - 1 - j].getImag() / 2;

                x2 = carr[i][csize2 - j].getReal() / 2;
                y2 = carr[i][csize2 - j].getImag() / 2;

                m_scratcha[2 * i      ] = complex<T>(x1 + x2, y1 - y2);
                m_scratcha[(2 * i) + 1] = complex<T>(y1 + y2, x2 - x1);

                m_scratchb[2 * i      ] = complex<T>(xb1 + xb2, yb1 - yb2);
                m_scratchb[(2 * i) + 1] = complex<T>(yb1 + yb2, xb2 - xb1);
            }

            for(i *= 2; i < m_vec_size1; ++i)
            {
                m_scratcha[i] = m_scratchb[i] = complex<T>(0., 0.);  /*!< Zero pad */
            }

            m_fft1.ForwardDecFreq(m_vec_size1, m_scratchb);
            m_fft1.ForwardDecFreq(m_vec_size1, m_scratcha);

            carr[0][j         ] = m_scratcha[0];
            carr[0][csize2 - j] = m_scratcha[0].conj();

            carr[0][j + 1         ] = m_scratchb[0];
            carr[0][csize2 - 1 - j] = m_scratchb[0].conj();

            for(i = 1; i < csize1; ++i)
            {
                carr[i][j] = m_scratcha[i];
                carr[i][j + 1] = m_scratchb[i];
                carr[i][csize2 - 1 - j] = m_scratchb[m_vec_size1 - i].conj();
                carr[i][csize2 - j] = m_scratcha[m_vec_size1 - i].conj();
            }
        }

        //! There should be 1 column left over.
        if(j < csize2 / 2)
        {
            for(i = 0; i < (rsize1 + 1) / 2; ++i)
            {
                x1 = carr[i][j].getReal() / 2;
                y1 = carr[i][j].getImag() / 2;

                x2 = carr[i][csize2 - j].getReal() / 2;
                y2 = carr[i][csize2 - j].getImag() / 2;

                m_scratcha[2 * i      ] = complex<T>(x1 + x2, y1 - y2);
                m_scratcha[(2 * i) + 1] = complex<T>(y1 + y2, x2 - x1);
            }

            for(i *= 2; i < m_vec_size1; ++i)
            {
                m_scratcha[i] = complex<T>(0., 0.); /*!< Zero pad */
            }

            m_fft1.ForwardDecFreq(m_vec_size1, m_scratcha);
            carr[0][j] = m_scratcha[0];
            carr[0][csize2-j] = m_scratcha[0].conj();

            for(i = 1; i < csize1; ++i)
            {
                carr[i][j] = m_scratcha[i];
                carr[i][csize2 - j] = m_scratcha[m_vec_size1 - i].conj();
            }
        }
    }

    /*!
     * \brief InverseRC - Inverse analogue of Forward FFT function ForwardRC. It is natural
     *                    to pair ForwardCR with InverseRC. The routine choice is made automatically
     *                    (based on speed estimates)by the routines ::Forward &
     *                    ::Inverse (below in the "public" access block).
     * \param csize1 - dimension of input array
     * \param csize2 - dimension of input array
     * \param carr - input 2 dimension array with complex coeff
     * \param rsize1 - dimension of output array
     * \param rsize2 - dimension of output array
     * \param rarr - output 2 dimension array
     */
    void InverseRC(int csize1,
                   int csize2,
                   const complex<T>* const* carr,
                   int rsize1,
                   int rsize2,
                   double** rarr)
    {
        SetupInverse(2 * (csize1 - 1), csize2); /*!< Safety */
        if(m_vec_size1 == 0 || m_vec_size2 == 0) return; /*!< Nothing to do */

        int i, j;
        double x1, y1, x2, y2;
        double xb1, yb1, xb2, yb2;

        //! Do row inverse FFT's. Handle the first & csize1'th row specially.
        //! These rows are the DFT's of real sequences, so they each satisfy the conjugate symmetry condition
        m_work_arr[0][0] = complex<T>(carr[0][0].getReal(), carr[csize1 - 1][0].getReal());

        for(j=1;j<csize2/2;j++)
        {
            x1 = carr[0][j].getReal();
            y1 = carr[0][j].getImag();

            x2 = carr[csize1 - 1][j].getReal();
            y2 = carr[csize1 - 1][j].getImag();

            m_work_arr[0][j         ] = complex<T>(x1 - y2, x2 + y1);
            m_work_arr[0][csize2 - j] = complex<T>(x1 + y2, x2 - y1);
        }

        m_work_arr[0][csize2 / 2] = complex<T>(carr[0][csize2 / 2].getReal(),
                                               carr[csize1 - 1][csize2 / 2].getReal());
        m_fft2.InverseDecTime(csize2, m_work_arr[0], 1.);

        //! iFFT the remaining rows
        for(i = 1; i < csize1 - 1; ++i)
        {
            for(j = 0; j < csize2; ++j)
            {
                m_work_arr[i][j] = carr[i][j];
            }
            m_fft2.InverseDecTime(csize2, m_work_arr[i], 1.);
        }

        // Now do iFFT's on columns.  These are conj. symmetric, so we
        // process them 2 at a time.  Also, recall the 1st row of workarr
        // contains the iFFT's of the 1st and csize1'th row of the given carr.
        for(j = 0; j + 3 < rsize2; j += 4)
        {
            m_scratcha[0] = complex<T>(m_work_arr[0][j].real(), m_work_arr[0][j + 1].real());
            m_scratcha[csize1 - 1] = complex<T>(m_work_arr[0][j].imag(), m_work_arr[0][j + 1].imag());
            m_scratchb[0] = complex<T>(m_work_arr[0][j + 2].real(), m_work_arr[0][j + 3].real());
            m_scratchb[csize1 - 1] = complex<T>(m_work_arr[0][j + 2].imag(), m_work_arr[0][j + 3].imag());

            for(i = 1; i < csize1 - 1; ++i)
            {
                x1 = m_work_arr[i][j].real();
                y1 = m_work_arr[i][j].imag();

                x2 = m_work_arr[i][j + 1].real();
                y2 = m_work_arr[i][j + 1].imag();

                xb1 = m_work_arr[i][j + 2].real();
                yb1 = m_work_arr[i][j + 2].imag();

                xb2 = m_work_arr[i][j + 3].real();
                yb2 = m_work_arr[i][j + 3].imag();

                m_scratcha[i]               = complex<T>(x1 - y2, x2 + y1);
                m_scratcha[m_vec_size1 - i] = complex<T>(x1 + y2, x2 - y1);
                m_scratchb[i]               = complex<T>(xb1 - yb2, xb2 + yb1);
                m_scratchb[m_vec_size1 - i] = complex<T>(xb1 + yb2, xb2 - yb1);
            }

            m_fft1.InverseDecTime(m_vec_size1, m_scratchb, static_cast<T>(m_vec_size1 * m_vec_size2));
            m_fft1.InverseDecTime(m_vec_size1, m_scratcha, static_cast<T>(m_vec_size1 * m_vec_size2));

            for(i = 0; i < rsize1; ++i)
            {
                rarr[i][j]     = m_scratcha[i].real();
                rarr[i][j + 1] = m_scratcha[i].imag();
                rarr[i][j + 2] = m_scratchb[i].real();
                rarr[i][j + 3] = m_scratchb[i].imag();
            }
        }

        // Remaining columns if rsize2 is not divisible by 4.  OTOH, csize2
        // *is* divisible by 2, so we can assume workarr[i][j+1] exists.
        for(; j < rsize2; j += 2)
        {
            m_scratcha[0] = complex<T>(m_work_arr[0][j].real(), m_work_arr[0][j + 1].real());
            m_scratcha[csize1 - 1] = complex<T>(m_work_arr[0][j].imag(), m_work_arr[0][j + 1].imag());

            for(i = 1; i < csize1 - 1; ++i)
            {
                x1 = m_work_arr[i][j].real();
                y1 = m_work_arr[i][j].imag();

                x2 = m_work_arr[i][j + 1].real();
                y2 = m_work_arr[i][j + 1].imag();

                m_scratcha[i]               = complex<T>(x1 - y2, x2 + y1);
                m_scratcha[m_vec_size1 - i] = complex<T>(x1 + y2, x2 - y1);
            }

            m_fft1.InverseDecTime(m_vec_size1, m_scratcha, static_cast<T>(m_vec_size1 * m_vec_size2));

            for(i=0;i<rsize1;i++)
            {
                rarr[i][j] = m_scratcha[i].real();
                if(j + 1 < rsize2) rarr[i][j + 1] = m_scratcha[i].imag();
            }
        }
    }

    /*!
     * \brief InverseCR - Inverse analogue of Forward FFT function ForwardCR. It is natural
     *                    to pair ForwardRC with InverseCR. The routine choice is made automatically
     *                    (based on speed estimates)by the routines ::Forward &
     *                    ::Inverse (below in the "public" access block).
     * \param csize1 - dimension of input array
     * \param csize2 - dimension of input array
     * \param carr - input 2 dimension array with complex coeff
     * \param rsize1 - dimension of output array
     * \param rsize2 - dimension of output array
     * \param rarr - output 2 dimension array
     */
    void InverseCR(int csize1,
                   int csize2,
                   const complex<T>* const* carr,
                   int rsize1,
                   int rsize2,
                   double** rarr)
    {
        SetupInverse(2*(csize1-1),csize2); /*!< Safety */
        if(m_vec_size1 == 0 || m_vec_size2 == 0) return; /*!< Nothing to do */

        int i, j;
        double x1, y1, x2, y2;
        double xb1, yb1, xb2, yb2;

        //! Column iFFT's Handle the first & csize2 / 2'th column specially.
        //! These cols are the DFT's of real sequences, so they each satisfy the conjugate
        //! symmetry condition
        m_scratcha[0] = complex<T>(carr[0][0].real(), carr[0][csize2 / 2].real());

        for(i = 1; i < csize1 - 1; ++i)
        {
            x1 = carr[i][0].real();
            y1 = carr[i][0].imag();

            x2 = carr[i][csize2 / 2].real();
            y2 = carr[i][csize2 / 2].imag();

            m_scratcha[i]            = complex<T>(x1 - y2, x2 + y1);
            m_scratcha[m_vec_size1 - i] = complex<T>(x1 + y2, x2 - y1);
        }

        m_scratcha[csize1 - 1] = complex<T>(carr[csize1 - 1][0].real(),
                                            carr[csize1 - 1][csize2 / 2].real());
        m_fft1.InverseDecTime(m_vec_size1, m_scratcha, 1);

        for(i = 0; i < m_vec_size1; i += 2)
        { //! ASSUMES vecsize1 is even
          //! See packing note below.
          m_work_arr[i / 2][0]          = complex<T>(m_scratcha[i].real(), m_scratcha[i+1].real());
          m_work_arr[i / 2][csize2 / 2] = complex<T>(m_scratcha[i].imag(), m_scratcha[i+1].imag());
        }

        //! Do remaining column iFFT's, two at a time for better memory access locality.
        for(j = 1; j + 1 < csize2 / 2; j += 2)
        {
            m_scratcha[0] = carr[0][j];
            m_scratchb[0] = carr[0][j + 1];

            for(i = 1; i < csize1 - 1; ++i)
            {
                m_scratcha[i] = carr[i][j];
                m_scratchb[i] = carr[i][j + 1];
                m_scratchb[m_vec_size1 - i] = carr[i][csize2 - 1 - j].conj();
                m_scratcha[m_vec_size1 - i] = carr[i][csize2 - j].conj();
            }
            m_scratcha[csize1 - 1] = carr[csize1 - 1][j];
            m_scratchb[csize1 - 1] = carr[csize1 - 1][j + 1];
            m_fft1.InverseDecTime(m_vec_size1, m_scratchb, 1.);
            m_fft1.InverseDecTime(m_vec_size1, m_scratcha, 1.);

            //! Pack into m_work_arr. Rows will be conjugate symmetric, so we
            //! can pack two rows into 1 via r[k]+i.r[k+1] -> m_work_arr[k/2].
            for(i = 0; i < rsize1; i += 2)
            {
                // CAREFUL! The above 'rsize1' bound may depend on how the
                // iFFT's are calculated in the 'Row iFFT's' code section,
                // and how 'i' is initialized.
                x1 = m_scratcha[i].real();
                y1 = m_scratcha[i].imag();

                x2 = m_scratcha[i + 1].real();
                y2 = m_scratcha[i + 1].imag();

                xb1 = m_scratchb[i].real();
                yb1 = m_scratchb[i].imag();

                xb2 = m_scratchb[i + 1].real();
                yb2 = m_scratchb[i + 1].imag();

                m_work_arr[i / 2][j]              = complex<T>(x1 - y2, x2 + y1);
                m_work_arr[i / 2][j + 1]          = complex<T>(xb1 - yb2, xb2 + yb1);
                m_work_arr[i / 2][csize2 - j - 1] = complex<T>(xb1 + yb2, xb2 - yb1);
                m_work_arr[i / 2][csize2 - j]     = complex<T>(x1 + y2, x2 - y1);
            }
        }

        //! There should be 1 column left over
        if((j = (csize2 / 2) - 1) % 2 == 1)
        {
            //! Column (csize2 / 2) - 1 *not* processed above
            m_scratcha[0] = carr[0][j];
            for(i = 1; i < csize1 - 1; ++i)
            {
                m_scratcha[i] = carr[i][j];
                m_scratcha[m_vec_size1 - i] = carr[i][csize2 - j].conj();
            }

            m_scratcha[csize1 - 1] = carr[csize1 - 1][j];
            m_fft1.InverseDecTime(m_vec_size1, m_scratcha, 1);

            for(i = 0; i < rsize1; i += 2)
            {
                // CAREFUL! The above 'rsize1' bound may depend on how the
                // iFFT's are calculated in the 'Row iFFT's' code section,
                // and how 'i' is initialized.
                x1 = m_scratcha[i].real();
                y1 = m_scratcha[i].imag();

                x2 = m_scratcha[i + 1].real();
                y2 = m_scratcha[i + 1].imag();

                m_work_arr[i / 2][j]          = complex<T>(x1 - y2, x2 + y1);
                m_work_arr[i / 2][csize2 - j] = complex<T>(x1 + y2, x2 - y1);
            }
        }

        //! Row iFFT's
        for(i = 0; i < rsize1; i += 2)
        {
            m_fft2.InverseDecTime(m_vec_size2, m_work_arr[i/2], static_cast<T>(m_vec_size1 * m_vec_size2));

            for(j = 0; j < rsize2; ++j)
            {
                rarr[i][j] = m_work_arr[i / 2][j].real();
            }

            if(i + 1 < rsize1)
            {
                for(j = 0; j < rsize2; ++j) rarr[i + 1][j] = m_work_arr[i / 2][j].imag();
            }
        }
    }

public:
    FFTR2D()
    {
        m_vec_size1 = m_vec_size2 = 0;
        m_log_size1 = m_log_size2 = 0;
        m_scratcha = nullptr;
        m_scratchb = nullptr;
        m_work_arr = nullptr;
    }

    ~FFTR2D()
    {
        ReleaseMemory();
    }

    // Forward 2D-FFT.
    // Computes FFT of rsize1 x rsize2 double** rarr, leaving result in
    // csize1 x csize2 MyComplex** carr.  rsize2 *must* be <= csize2, and
    // rsize1 *must* be <= 2*(csize1-1).  This routine returns only the
    // top half +1 of the transform.  The bottom half is given by
    //      carr[2*(csize1-1)-i][csize2-j]=conj(carr[i][j])
    // for i>=csize1, with the second indices interpreted 'mod csize2'.
    // The client can pass a full-sized top-half-filled MyComplex** carr
    // arr to the FFTReal2D member function "FillOut" to perform
    // this operation...
    //  NOTE: The Microsoft Visual C++ compiler apparently doesn't want to
    // automatically convert from (double**) to (const double* const*), so
    // you will need to put in an explicit cast for that compiler (and it
    // does no harm in any case).
    void Forward(int rsize1, int rsize2,
                 const double* const* rarr,
                 int csize1, int csize2, complex<T>** carr)
    {
        if(csize2<rsize2)
          PlainError(1,"Error in FFTRealD::Forward(int,int,double**,int,int,MyComplex**): "
               "csize2 (=%d) *must* be >= rsize2 (=%d)\n",csize2,rsize2);

        if(csize1<(rsize1/2)+1)
          PlainError(1,"Error in FFTRealD::Forward(int,int,double**,int,int,MyComplex**): "
               "csize1 (=%d) *must* be >= (rsize1/2)+1 (=%d)\n",
               csize1,(rsize1/2)+1);

        Setup(2*(csize1-1),csize2);
        if(vecsize1==0 || vecsize2==0) return; // Nothing to do

        // Determine which Forward routine to call (ForwardCR or ForwardRC)
        // (Use double arithmetic to protect against integer overflow.  If
        // the two times are very close, then the choice doesn't really
        // matter.)
        // 1) Estimated (proportional) time for ForwardCR
        double crtime=double(vecsize1*vecsize2)*double(logsize2)
                    + double(vecsize1*rsize2)*double(logsize1);
        // 2) Estimated (proportional) time for ForwardRC
        double rctime=double(rsize1*vecsize2)*double(logsize2)
                    + double(vecsize1*vecsize2)*double(logsize1);
        // Introduce empirical adjustment factor
        rctime*=CRRCspeedratio;

        if(crtime<=rctime) ForwardCR(rsize1,rsize2,rarr,csize1,csize2,carr);
        else               ForwardRC(rsize1,rsize2,rarr,csize1,csize2,carr);
    }

    // See previous note.  Here csize1 is the full array height, and the
    // top (csize1/2)+1 rows are assumed already filled.
    void FillOut(int csize1, int csize2, complex<T>** carr)
    {
        // This routine assumes carr is a top-half-filled DFT of
        // a real function, and fills in the bottom half using the
        // relation
        //      carr[csize1-i][csize2-j]=conj(carr[i][j])
        // for i>csize1/2, with the second indices interpreted 'mod csize2'.
        for(auto i = 1; i < csize1 / 2; ++i)
        {
            carr[csize1 - i][0] = carr[i][0]->conj();
            for(auto j = 1; j < csize2; ++j)
            {
                carr[csize1-i][j] = carr[i][csize2-j]->conj();
            }
        }
    }


    // Inverse 2D-FFT.  See notes to '::Forward', including the
    // const* note for MVC++.
    void Inverse(int csize1, int csize2,
                 const complex<T>* const* carr,
                 int rsize1, int rsize2, double** rarr)
    {
        if(csize2<rsize2)
          PlainError(1,"Error in FFTRealD::Inverse(int,int,double**,int,int,MyComplex**): "
               "csize2 (=%d) *must* be >= rsize2 (=%d)\n",csize2,rsize2);

        if(csize1<(rsize1/2)+1)
          PlainError(1,"Error in FFTRealD::Inverse(int,int,double**,int,int,MyComplex**): "
               "csize1 (=%d) *must* be >= (rsize1/2)+1 (=%d)\n",
               csize1,(rsize1/2)+1);

        SetupInverse(2*(csize1-1),csize2);
        if(vecsize1==0 || vecsize2==0) return; // Nothing to do

        // Determine which Inverse routine to call (InverseRC or InverseCR)
        // (Use double arithmetic to protect against integer overflow.  If
        // the two times are very close, then the choice doesn't really
        // matter.)
        // 1) Estimated (proportional) time for InverseRC (==ForwardCR)
        double irctime=double(vecsize1*vecsize2)*double(logsize2)
                    + double(vecsize1*rsize2)*double(logsize1);
        // 2) Estimated (proportional) time for InverseCR (==ForwardRC)
        double icrtime=double(rsize1*vecsize2)*double(logsize2)
                    + double(vecsize1*vecsize2)*double(logsize1);
        // Introduce empirical adjustment factor
        icrtime*=CRRCspeedratio;

        if(irctime<=icrtime) InverseRC(csize1,csize2,carr,rsize1,rsize2,rarr);
        else                 InverseCR(csize1,csize2,carr,rsize1,rsize2,rarr);
    }
};

} //namespace adf

#endif //FFT_H
