#include "test.h"


namespace test {

    void print_real(const std::vector<complex> &x) {
        for (const complex &val: x) {
            std::cout << val.real() << ", ";
        }
        std::cout << '\n';
    }

    void print_imag(const std::vector<complex> &x) {
        for (const complex &val: x) {
            std::cout << val.imag() << ", ";
        }
        std::cout << '\n';
    }

    void print_abs(const std::vector<complex> &x) {
        for (const complex &val: x) {
            std::cout << std::abs(val) << ", ";
        }
        std::cout << '\n';
    }

    std::vector<complex> dft(const std::vector<complex> &x) {
        int N = (int) x.size();
        std::vector<std::complex<double>> spectraComplex(N);

        for (int n = 0; n != N; ++n) {
            for (int k = 0; k != N; ++k) {
                spectraComplex[n] += x[k] * std::polar(1., -2 * M_PI * k * n / N);
            }
        }
        return spectraComplex;
    }


    void test_fft() {
        std::vector<complex> x2{1, 2, 3, 4, 5, 6, 7, 8};
        std::vector<complex> x3{1, 2, 3, 4, 5, 6, 7, 8, 9};
        std::vector<complex> x5(25);
        for (int i = 0; i < 25; ++i) {
            x5[i] = i + 1;
        }
        std::vector<complex> x_mix(30);
        for (int i = 0; i < 30; ++i) {
            x_mix[i] = i + 1;
        }
        std::vector<complex> out2 = FFT::fft(x2, false);
        std::vector<complex> out3 = FFT::fft(x3, false);
        std::vector<complex> out5 = FFT::fft(x5, false);
        std::vector<complex> out_mix = FFT::fft(x_mix, false);

        std::vector<complex> x2_after = FFT::fft(out2, true);
        std::vector<complex> x3_after = FFT::fft(out3, true);
        std::vector<complex> x5_after = FFT::fft(out5, true);
        std::vector<complex> x_mix_after = FFT::fft(out_mix, true);

        std::cout << "длина 2^3" << '\n';
        print_abs(out2);
        print_abs(x2_after);
        std::cout << "длина 3^2" << '\n';
        print_abs(out3);
        print_abs(x3_after);
        std::cout << "длина 5^2" << '\n';
        print_abs(out5);
        print_abs(x5_after);
        std::cout << "длина 2 * 3 * 5" << '\n';
        print_abs(out_mix);
        print_abs(x_mix_after);

    }

}