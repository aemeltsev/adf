# adf - :alembic: c++ analog and digital filter library (In Progress)

The repository contains the adf library code. Adf is an open-source solution of C++ code (compliant to C++11) designed for synthesis and design filters.

# About adf

Adf is an open-source C++ classes collection library designed to implementation of analog active and digital iir and fir filters.
Inside adf includes estimate normalized transfer functions for the Butterworth, Chebyshev, inverse Chebyshev, and elliptic approximation case. Conversion of the normalized low pass filter to an denormalized low pass, high pass, band pass, or band stop filter, and addition, the calculation of the frequency response for analog filter. For design digital iir (recursive) filter using method of a bilinear z-transform. With further calculation the frequency response. Also in library contained custom fast fourier transformation module.

# Quickstart
todo

# Building adf

To build adf you need using CMake.

# Codemap

Adf contains the following C++ modules components:

* advmath Advanced mathematical methods.
<br /> Contain hyperbolic, elliptic and other special methods which all other adf code depends on.

* anfilter Calculate of the analog active filters

* base Internal basic methods and definitions. Not used outside.

* complex Basic implementation class of the complex number.

* fft Fast Fourier transformation.
 
* fresp Calculate the frequency response of a filter using starting and stopping frequency.

* iirfilter Calculate of the digital iir-filter.

* polynomial Custom representation of the polynomial class.

# License

The adf c++ library is licensed under the terms of the Apache license. See LICENSE for more information.

# References

For more information about adf and methods for filter design, see:

Practical Analog and Digital Filter Design, Artech House, *L. Thede*

Analog Electronic Filters. Theory, Design and Synthesis, 2012th ed. Springer, *H.G. Dimopoulos*

Passive, Active, and Digital Filters, 2nd ed. CRC Press, *Wai-Kai Chen*

Modern Filter theory and design, John Wiley & Sons, *G. C. Temes*

Digital Signal Processing. A Practical Approach, 2nd ed. Prentice Hall, *E. C. Ifeachor*
