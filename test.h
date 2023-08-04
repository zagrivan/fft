#ifndef FFT_TEST_H
#define FFT_TEST_H

#include "fft.h"

namespace test {
    void print_real(const std::vector<complex>& x);
    void print_imag(const std::vector<complex>& x);
    void print_abs(const std::vector<complex>& x);
    std::vector<complex> dft(const std::vector<complex> &x);
}


#endif //FFT_TEST_H
