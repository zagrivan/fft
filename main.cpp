#include <iostream>
#include <complex>
#include <vector>

#include "fft.h"
#include "test.h"


int main(int argc, char *argv[]) {

    std::vector<complex> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,27,28};
    auto out = FFT::fft(v, false);
    test::print_abs(out);
    auto out1 = FFT::fft(out, true);
    test::print_abs(out);
    return 0;
}
