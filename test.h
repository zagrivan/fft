#ifndef FFT_TEST_H
#define FFT_TEST_H

#include "fft.h"

namespace test {
    void print_real(const std::vector<complex>& x);
    void print_imag(const std::vector<complex>& x);
    void print_abs(const std::vector<complex>& x);
    std::vector<complex> dft(const std::vector<complex> &x);
    void test_fft();
    std::vector<complex> get_cos(int n, double f, double Ts);
    std::vector<complex> get_sin_t_t(int n, double f, double Ts);
}


#endif //FFT_TEST_H
