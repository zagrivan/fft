#ifndef FFT_FFT_H
#define FFT_FFT_H

#include <complex>
#include <vector>
#include <iostream>


typedef std::complex<double> complex;

constexpr complex J(0, 1);
// для radix3
const double sinPi3 = std::sin(M_PI / 3.);
// для radix5
const complex w51(cos(-2. * M_PI / 5), sin(-2. * M_PI / 5));  // exp(-2πi*1/5)
const complex w52(cos(-4. * M_PI / 5), sin(-4. * M_PI / 5)); // exp(-2πi*2/5)
const complex w53(cos(-6. * M_PI / 5), sin(-6. * M_PI / 5)); // exp(-2πi*3/5)
const complex w54(cos(-8. * M_PI / 5), sin(-8. * M_PI / 5));  // exp(-2πi*4/5)


class FFT {
public:
    static void fft(std::vector<complex> &in, bool inverse=false);
private:
    FFT()= default;
    ~FFT()= default;

    static int validate(int n);

    static void radix2(std::vector<complex> &x, double direction);

    static void radix3(std::vector<complex> &x, double direction);

    static void radix5(std::vector<complex> &x, double direction);

    static void radix_mix(std::vector<complex> &x, double direction);

};

#endif //FFT_FFT_H
