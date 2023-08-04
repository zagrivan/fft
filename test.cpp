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

}