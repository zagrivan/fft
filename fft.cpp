#include <vector>
#include <complex>

#include "fft.h"

int FFT::validate(int n) {
    while (n % 2 == 0) n /= 2;
    while (n % 3 == 0) n /= 3;
    while (n % 5 == 0) n /= 5;
    if (n == 1)
        return 30; // если число n состоит из нескольких простых чисел
    return -1;  // если одно из простых чисел больше 5
}

void FFT::radix2(std::vector<complex> &x, double direction) {
    int N = (int) x.size();

    std::vector<complex> p0(N / 2), p1(N / 2);

    for (int i = 0; i < N / 2; ++i) {
        p0[i] = x[2 * i];
        p1[i] = x[2 * i + 1];
    }

    radix_mix(p0, direction);
    radix_mix(p1, direction);

    complex w = 1;
    complex wn = std::polar(1., direction * 2 * M_PI / N);

    for (int i = 0; i < N / 2; ++i) {
        x[i] = p0[i] + w * p1[i];
        x[i + N / 2] = p0[i] - w * p1[i];
        w *= wn;
    }
}

void FFT::radix3(std::vector<complex> &x, double direction) {
    int N = (int) x.size();

    std::vector<complex> p0(N / 3), p1(N / 3), p2(N / 3);

    for (int i = 0; i < N / 3; ++i) {
        p0[i] = x[3 * i];
        p1[i] = x[3 * i + 1];
        p2[i] = x[3 * i + 2];
    }

    radix_mix(p0, direction);
    radix_mix(p1, direction);
    radix_mix(p2, direction);

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

void FFT::radix5(std::vector<complex> &x, double direction) {
    int N = (int) x.size();

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

    complex w = 1;
    complex wn = std::polar(1., direction * 2 * M_PI / N);

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
        } else {              // обратное
            x[i + N / 5] = t0 + t1 * w54 + t2 * w53 + t3 * w52 + t4 * w51;
            x[i + 2 * N / 5] = t0 + t1 * w53 + t2 * w51 + t3 * w54 + t4 * w52;
            x[i + 3 * N / 5] = t0 + t1 * w52 + t2 * w54 + t3 * w51 + t4 * w53;
            x[i + 4 * N / 5] = t0 + t1 * w51 + t2 * w52 + t3 * w53 + t4 * w54;
        }
        w *= wn;
    }
}

// длина последовательности кратна нескольким простым числам(2, 3, 5), fft вычисляется рекурсивно
void FFT::radix_mix(std::vector<complex> &x, double direction) {
    int N = (int) x.size();
    if (N == 1) return;

    if (N % 2 == 0) {
        radix2(x, direction);
    } else if (N % 3 == 0) {
        radix3(x, direction);
    } else if (N % 5 == 0) {
        radix5(x, direction);
    }
}

void FFT::fft(std::vector<complex> &in, bool inverse) {
    int N = (int) in.size();
    if (N < 2) return;
    double direction = inverse ? 1 : -1; // 1 - обратное, -1 - прямое

    if (validate(N) == -1) {
        std::cerr << "The length of the sequence must be a multiple of only 2, 3 and 5\n";
        throw std::exception();
    }

    FFT::radix_mix(in, direction); /// вычисление fft

    if (inverse) {                   // для обратного fft умножение на 1/N
        for (complex &val: in) {
            val /= N;
        }
    }
}