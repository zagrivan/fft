#include <vector>
#include <complex>

#include "fft.h"


int FFT::radix(int n) {
    if ((n & (n - 1)) == 0) {
        return 2;
    }
    for (int i = 3; i <= 5; i += 2) {
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
    if (n == 1) return 30;
    return -1;
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

std::vector<complex> FFT::fft_2_iterative(const std::vector<complex> &xt, double direction) {
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
        complex wn = std::polar(1., direction * 2 * M_PI / m);
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

void FFT::fft_radix3_rec(std::vector<complex> &x, double direction) {
    int N = (int) x.size();
    if (N == 1) return;

    std::vector<complex> p0(N / 3), p1(N / 3), p2(N / 3);

    for (int i = 0; i < N / 3; ++i) {
        p0[i] = x[3 * i];
        p1[i] = x[3 * i + 1];
        p2[i] = x[3 * i + 2];
    }

    fft_radix3_rec(p0, direction);
    fft_radix3_rec(p1, direction);
    fft_radix3_rec(p2, direction);

    complex w = 1;
    complex wn = std::polar(1., direction * 2 * M_PI / N);

    for (int i = 0; i < N / 3; ++i) {
        complex t0 = w * p1[i] + w * w * p2[i];
        complex t1 = p0[i] - t0 / 2.0;
        complex t2 = J * sinPi3 * (w * p1[i] - w * w * p2[i]) * direction;

        x[i] = p0[i] + t0;
        x[i + N / 3] = t1 + t2;
        x[i + 2 * N / 3] = t1 - t2;

        w *= wn;
    }
}


void FFT::fft_radix5_rec(std::vector<complex> &x, double direction) {
    int N = (int) x.size();
    if (N == 1) return;

    std::vector<complex> p0(N / 5), p1(N / 5), p2(N / 5), p3(N / 5), p4(N / 5);

    for (int i = 0; i < N / 5; ++i) {
        p0[i] = x[5 * i];
        p1[i] = x[5 * i + 1];
        p2[i] = x[5 * i + 2];
        p3[i] = x[5 * i + 3];
        p4[i] = x[5 * i + 4];
    }

    fft_radix5_rec(p0, direction);
    fft_radix5_rec(p1, direction);
    fft_radix5_rec(p2, direction);
    fft_radix5_rec(p3, direction);
    fft_radix5_rec(p4, direction);

    complex w = 1;
    complex wn = std::polar(1., direction * 2 * M_PI / N);

    for (int i = 0; i < N; ++i) {
        int index = i % (N / 5);
        x[i] = p0[index] + p1[index] * w + p2[index] * std::pow(w, 2) + p3[index] * std::pow(w, 3) + p4[index] * std::pow(w, 4);
        w *= wn;
    }
}


void FFT::fft_mix_radix(std::vector<complex> &x, double direction) {
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

        fft_mix_radix(p0, direction);
        fft_mix_radix(p1, direction);

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

        fft_mix_radix(p0, direction);
        fft_mix_radix(p1, direction);
        fft_mix_radix(p2, direction);

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

        fft_mix_radix(p0, direction);
        fft_mix_radix(p1, direction);
        fft_mix_radix(p2, direction);
        fft_mix_radix(p3, direction);
        fft_mix_radix(p4, direction);

        for (int i = 0; i < N; ++i) {
            int index = i % (N / 5);
            x[i] = p0[index] + p1[index] * w + p2[index] * w * w + p3[index] * std::pow(w, 3) + p4[index] * std::pow(w, 4);
            w *= wn;
        }
    }
}


std::vector<complex> FFT::fft(const std::vector<complex> &in, bool inverse) {
    int N = (int) in.size();
    if (N == 0) return {};
    int rdx = FFT::radix(N);
    double direction = inverse ? 1 : -1;

    if (rdx == -1) {
        std::cerr << "The length of the sequence must be a multiple of only 2, 3 and 5\n";
        throw std::exception();
    }

    std::vector<complex> out;

    if (rdx == 2) {
        out = FFT::fft_2_iterative(in, direction);
    } else {
        out = in;
        if (rdx == 3) {
            FFT::fft_radix3_rec(out, direction);
        } else if (rdx == 5) {
            FFT::fft_radix5_rec(out, direction);
        } else if (rdx == 30) {
            FFT::fft_mix_radix(out, direction);
        }
    }

    if (inverse) {
        for (complex &val: out) {
            val /= N;
        }
    }

    return out;
}