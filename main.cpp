#include <complex>
#include <vector>

#include "fft.h"
#include "test.h"


int main() {
    // длина последовательности 2^9 = 512
    int n = (int) std::pow(2, 9);
    // генерируем n отсчетов, представляющих сумму 4 cos с минимальной частотой f=0.3
    std::vector<complex> seq = test::get_cos(n, 0.3, 10.0 / n);

    test::print_complex(seq); // перед преобразованием

    FFT::fft(seq, false); // forward fft

    test::print_complex(seq); // после прямого

    FFT::fft(seq, true); // backward fft

    test::print_complex(seq); // после обратного

    return 0;
}
