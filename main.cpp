#include <iostream>
#include <complex>
#include <vector>

#include "fft.h"
#include "test.h"


int main(int argc, char *argv[]) {
    // длина последовательности 2^9 = 512
    int n = (int) std::pow(2, 9);
    // генерируем n отсчетов, представляющих сумму 4 cos с минимальной частотой f=0.3
    std::vector<complex> in = test::get_cos(n, 0.3, 10.0 / n);

    std::vector<complex> after_fft = FFT::fft(in, false); // forward fft

    std::vector<complex> after_inverse_fft = FFT::fft(after_fft, true); // backward fft

    test::print_real(in);
    test::print_abs(after_fft);  // вывод комплексного числа по модулю( sqrt( real^2 + imag^2 ) )
    test::print_real(after_inverse_fft);

    return 0;
}
