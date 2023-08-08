#ifndef FFT_FFT_H
#define FFT_FFT_H


#include <complex>
#include <vector>
#include <iostream>


typedef std::complex<double> complex;

const complex J(0, 1);
const double sinPi3 = std::sin(M_PI / 3.); // для radix3
// для radix5
const complex w51(cos(-2. * M_PI / 5), sin(-2. * M_PI / 5));  // exp(-2πi*1/5)
const complex w52(cos(-4. * M_PI / 5), sin(-4. * M_PI / 5)); // exp(-2πi*2/5)
const complex w53(cos(-6. * M_PI / 5), sin(-6. * M_PI / 5)); // exp(-2πi*3/5)
const complex w54(cos(-8. * M_PI / 5), sin(-8. * M_PI / 5));  // exp(-2πi*4/5)



class FFT {
public:
    static std::vector<complex> fft(const std::vector<complex> &in, bool inverse);
private:
    FFT()= default;

    static int radix(int n);

    static void add_zeros(std::vector<complex> &p);

    static void reorder_base_r(std::vector<complex> &x, int radix, int n);

    static unsigned int bitReverse(unsigned int x, int log2n);

    static std::vector<complex> radix2(const std::vector<complex> &xt, double direction);

    static std::vector<complex> radix3(std::vector<complex> x, double direction);

    static std::vector<complex> radix5(std::vector<complex> x, double direction);

    static void radix_mix(std::vector<complex> &x, double direction);

};

#endif //FFT_FFT_H
