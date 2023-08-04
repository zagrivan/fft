#ifndef FFT_FFT_H
#define FFT_FFT_H


#include <complex>
#include <vector>
#include <iostream>


typedef std::complex<double> complex;

const complex J(0, 1);
const double sinPi3 = std::sin(M_PI / 3.);


class FFT {
public:
    static std::vector<complex> fft(const std::vector<complex> &in, bool inverse);
private:
    FFT()= default;

    static int radix(int n);

    static void add_zeros(std::vector<complex> &p);

    static unsigned int bitReverse(unsigned int x, int log2n);

    static void fft_mix_radix(std::vector<complex> &x, double direction);

    static std::vector<complex> fft_2_iterative(const std::vector<complex> &xt, double direction);

    static void fft_radix3_rec(std::vector<complex> &x, double direction);

    static void fft_radix5_rec(std::vector<complex> &x, double direction);

};

#endif //FFT_FFT_H
