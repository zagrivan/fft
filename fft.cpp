#include <vector>
#include <complex>

#include "fft.h"


int FFT::radix(int n) {
    if ((n & (n - 1)) == 0) {
        return 2;
    }
    for (int i = 3; i <= 5; i += 2) { // проверка для 3 и 5
        int n1 = n;
        while (n1 % i == 0) {
            n1 /= i;
        }
        if (n1 == 1)
            return i;
    }
    while (n % 2 == 0) n /= 2;
    while (n % 3 == 0) n /= 3;
    while (n % 5 == 0) n /= 5;
    if (n == 1)
        return 30; // если число n состоит из нескольких простых чисел
    return -1;  // если одно из простых чисел больше 5
}

void FFT::add_zeros(std::vector<complex> &p) {
    while ((p.size() & (p.size() - 1)) != 0) {
        p.emplace_back(0);
    }
}

unsigned int FFT::bitReverse(unsigned int x, int log2n) {
    int n = 0;
    for (int i = 0; i < log2n; i++) {
        n <<= 1;
        n |= x & 1;
        x >>= 1;
    }
    return n;
}

// для итеративных версий radix 3 и 5, перестанавливает входные значения для fft
void FFT::reorder_base_r(std::vector<complex> &x, int radix, int n) {
    for (int i = 0, j = 0; i < n - 1; ++i) {
        if (i < j) {
            std::swap(x[i], x[j]);
        }
        int k = (radix - 1) * n / radix;
        while (k <= j) {
            j -= k;
            k /= radix;
        }
        j += k / (radix - 1);
    }
}


std::vector<complex> FFT::radix2(const std::vector<complex> &xt, double direction) {
    int n = (int) xt.size();
    int log2n = (int) std::log2(n);
    std::vector<complex> Xf(n);
    for (unsigned int i = 0; i < n; ++i) {
        Xf[bitReverse(i, log2n)] = xt[i];
    }
    for (int s = 1; s <= log2n; ++s) {
        int m = 1 << s; // 2 ** s
        int m2 = m >> 1; // 2 ** s / 2
        complex w(1, 0);
        complex wn = std::polar(1., direction * 2 * M_PI / m); // exp(-1 * 2*PI/m) - прямое, exp(1 * 2*PI/m) - обратное
        for (int j = 0; j < m2; ++j) {
            for (int k = j; k < n; k += m) {
                complex a = Xf[k];
                complex b = w * Xf[k + m2];
                Xf[k] = a + b;
                Xf[k + m2] = a - b;
            }
            w *= wn;
        }
    }
    return Xf;
}

std::vector<complex> FFT::radix3(std::vector<complex> x, double direction) {
    const int r = 3;
    int n = (int) x.size();

    reorder_base_r(x, r, n);

    for (int m = r; m <= n; m *= r) {
        complex w(1, 0);
        complex wn = std::polar(1., direction * 2 * M_PI / m);
        for (int i = 0; i < m / r; ++i) {
            for (int j = 0; j < n; j += m) {
                complex t0 = x[i + j + m / r] * w + x[i + j + 2 * m / r] * w * w;
                complex t1 = x[i + j] - t0 / 2.;
                complex t2 = J * sinPi3 * (x[i + j + m / r] * w - x[i + j + 2 * m / r] * w * w) * direction;

                x[i + j] += t0;
                x[i + j + m / r] = t1 + t2;
                x[i + j + 2 * m / r] = t1 - t2;
            }
            w *= wn;
        }
    }

    return x;
}

std::vector<complex> FFT::radix5(std::vector<complex> x, double direction) {
    const int r = 5;
    int n = (int) x.size();

    reorder_base_r(x, r, n);

    for (int m = r; m <= n; m *= r) {
        complex w(1, 0);
        complex wn = std::polar(1., direction * 2 * M_PI / m);
        for (int i = 0; i < m / r; ++i) {
            for (int j = 0; j < n; j += m) {
                complex t0 = x[i + j];
                complex t1 = x[i + j + m / r] * w;
                complex t2 = x[i + j + 2 * m / r] * w * w;
                complex t3 = x[i + j + 3 * m / r] * w * w * w;
                complex t4 = x[i + j + 4 * m / r] * w * w * w * w;

                x[i + j] += t1 + t2 + t3 + t4;
                if (direction < 0) {  // прямое
                    x[i + j + m / r] = t0 + t1 * w51 + t2 * w52 + t3 * w53 + t4 * w54;
                    x[i + j + 2 * m / r] = t0 + t1 * w52 + t2 * w54 + t3 * w51 + t4 * w53;
                    x[i + j + 3 * m / r] = t0 + t1 * w53 + t2 * w51 + t3 * w54 + t4 * w52;
                    x[i + j + 4 * m / r] = t0 + t1 * w54 + t2 * w53 + t3 * w52 + t4 * w51;
                } else {    // обратное
                    x[i + j + m / r] = t0 + t1 * w54 + t2 * w53 + t3 * w52 + t4 * w51;
                    x[i + j + 2 * m / r] = t0 + t1 * w53 + t2 * w51 + t3 * w54 + t4 * w52;
                    x[i + j + 3 * m / r] = t0 + t1 * w52 + t2 * w54 + t3 * w51 + t4 * w53;
                    x[i + j + 4 * m / r] = t0 + t1 * w51 + t2 * w52 + t3 * w53 + t4 * w54;
                }
            }
            w *= wn;
        }
    }
    return x;
}

// для случая, когда длина последовательности кратна нескольким простым числам(2, 3, 5), fft вычисляется рекурсивно
void FFT::radix_mix(std::vector<complex> &x, double direction) {
    int N = (int) x.size();
    if (N == 1) return;

    complex w = 1;
    complex wn = std::polar(1., direction * 2 * M_PI / N);

    if (N % 2 == 0) {
        std::vector<complex> p0(N / 2), p1(N / 2);

        for (int i = 0; i < N / 2; i++) {
            p0[i] = x[2 * i];
            p1[i] = x[2 * i + 1];
        }

        radix_mix(p0, direction);
        radix_mix(p1, direction);

        for (int i = 0; i < N / 2; ++i) {
            x[i] = p0[i] + w * p1[i];
            x[i + N / 2] = p0[i] - w * p1[i];
            w *= wn;
        }
    } else if (N % 3 == 0) {
        std::vector<complex> p0(N / 3), p1(N / 3), p2(N / 3);

        for (int i = 0; i < N / 3; ++i) {
            p0[i] = x[3 * i];
            p1[i] = x[3 * i + 1];
            p2[i] = x[3 * i + 2];
        }

        radix_mix(p0, direction);
        radix_mix(p1, direction);
        radix_mix(p2, direction);

        for (int i = 0; i < N / 3; ++i) {
            complex t0 = w * p1[i] + w * w * p2[i];
            complex t1 = p0[i] - t0 / 2.0;
            complex t2 = J * sinPi3 * (w * p1[i] - w * w * p2[i]) * direction;

            x[i] = p0[i] + t0;
            x[i + N / 3] = t1 + t2;
            x[i + 2 * N / 3] = t1 - t2;

            w *= wn;
        }
    } else if (N % 5 == 0) {
        std::vector<complex> p0(N / 5), p1(N / 5), p2(N / 5), p3(N / 5), p4(N / 5);

        for (int i = 0; i < N / 5; ++i) {
            p0[i] = x[5 * i];
            p1[i] = x[5 * i + 1];
            p2[i] = x[5 * i + 2];
            p3[i] = x[5 * i + 3];
            p4[i] = x[5 * i + 4];
        }

        radix_mix(p0, direction);
        radix_mix(p1, direction);
        radix_mix(p2, direction);
        radix_mix(p3, direction);
        radix_mix(p4, direction);

        for (int i = 0; i < N / 5; ++i) {
            complex t0 = p0[i];
            complex t1 = p1[i] * w;
            complex t2 = p2[i] * w * w;
            complex t3 = p3[i] * w * w * w;
            complex t4 = p4[i] * w * w * w * w;

            x[i] = t0 + t1 + t2 + t3 + t4;
            if (direction < 0) {  // прямое
                x[i + N / 5] = t0 + t1 * w51 + t2 * w52 + t3 * w53 + t4 * w54;
                x[i + 2 * N / 5] = t0 + t1 * w52 + t2 * w54 + t3 * w51 + t4 * w53;
                x[i + 3 * N / 5] = t0 + t1 * w53 + t2 * w51 + t3 * w54 + t4 * w52;
                x[i + 4 * N / 5] = t0 + t1 * w54 + t2 * w53 + t3 * w52 + t4 * w51;
            } else {    // обратное
                x[i + N / 5] = t0 + t1 * w54 + t2 * w53 + t3 * w52 + t4 * w51;
                x[i + 2 * N / 5] = t0 + t1 * w53 + t2 * w51 + t3 * w54 + t4 * w52;
                x[i + 3 * N / 5] = t0 + t1 * w52 + t2 * w54 + t3 * w51 + t4 * w53;
                x[i + 4 * N / 5] = t0 + t1 * w51 + t2 * w52 + t3 * w53 + t4 * w54;
            }

            w *= wn;
        }
    }
}

std::vector<complex> FFT::fft(const std::vector<complex> &in, bool inverse) {
    int N = (int) in.size();
    if (N < 2) return in;
    int rdx = FFT::radix(N);
    double direction = inverse ? 1 : -1; // 1 - обратное, -1 - прямое

    if (rdx == -1) {
        std::cerr << "The length of the sequence must be a multiple of only 2, 3 and 5\n";
        throw std::exception();
    }

    std::vector<complex> out;

    if (rdx == 2) {
        out = FFT::radix2(in, direction);
    }
    else if (rdx == 3) {
        out = FFT::radix3(in, direction);
    }
    else if (rdx == 5) {
        out = FFT::radix5(in, direction);
    }
    else if (rdx == 30) {
        out = in;
        FFT::radix_mix(out, direction);
    }

    if (inverse) {
        for (complex &val: out) {
            val /= N;
        }
    }

    return out;
}