# fft
### Класс, реализующий быстрое преобразование фурье для длинны преобразования кратной 2, 3, 5.

В функцию `fft()`, определенную в классе `FFT`,
```c++
typedef std::complex<double> complex;
std::vector<complex> FFT::fft(const std::vector<complex> &in, bool inverse);
```
передается входная последовательность .
Для прямого преобразования `inverse = false`, для обратного `inverse = true`. Длина последовательности должна быть составным числом, кратным только 2, 3, 5.

В качестве примера работы в функции main файла `main.cpp` представлено прямое и обратное fft над последовательностью, представляющей сумму четырех косинусов.
```c++
// длина последовательности 2^9 = 512
int n = (int) std::pow(2, 9);
// генерируем n отсчетов, представляющих сумму 4 cos с минимальной частотой f=0.3
std::vector<complex> in = test::get_cos(n, 0.3, 10.0 / n);
std::vector<complex> after_fft = FFT::fft(in, false); // forward fft
std::vector<complex> after_inverse_fft = FFT::fft(after_fft, true); // backward fft
// вывод результата
test::print_real(in);
test::print_abs(after_fft);  // вывод комплексного числа по модулю( sqrt( real^2 + imag^2 ) )
test::print_real(after_inverse_fft);
```
Визуализация данных после обратного fft(равна последовательности до прямого fft)
![data](images/spectr.png "Визуализация данных после обратного fft")
Визуализация данных после прямого fft
![spectr](images/cos.png "Визуализация данных после прямого fft")
