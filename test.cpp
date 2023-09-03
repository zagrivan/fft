#include "test.h"


namespace test {

    void print_complex(const std::vector<complex> &x) {
        for (const complex& c : x) {
            std::cout << c << ", ";
        }
        std::cout << '\n';
    }

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

    std::vector<complex> get_cos(int n, double f, double Ts) {
        std::vector<complex> out(n);
        for (int i = 0; i < n; ++i) {
            out[i] = 2 * std::cos(2. * M_PI * f * i * Ts) + 0.9 * std::cos(2. * M_PI * f * 3 * i * Ts) + 0.6 * std::cos(2. * M_PI * f * 9 * i * Ts) + 0.4 * std::cos(2. * M_PI * f * 40 * i * Ts);
        }
        return out;
    }

    std::vector<complex> get_sin_t_t(int n, double f, double Ts) {
        std::vector<complex> out(n);
        for (int i = -n / 2, j = 0; i < n / 2; ++i, ++j) {
            if (i == 0) {
                out[j] = out[j - 1];
                continue;
            }
            out[j] = std::sin(2. * M_PI * f * i * Ts) / (i * Ts);
        }
        return out;
    }

}