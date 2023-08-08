# fft
### Класс, реализующий быстрое преобразование фурье для длинны преобразования кратной 2, 3, 5.

В функцию `fft()`, определенную в классе `FFT`,
```c++
typedef std::complex<double> complex;
std::vector<complex> FFT::fft(const std::vector<complex> &in, bool inverse);
```
передается входная последовательность .
Для прямого преобразования `inverse = false`, для обратного `inverse = true`. Длина последовательности должна быть составным числом, кратным только 2, 3, 5.

В качестве примера работы `FFT::fft()` функция `test::test_fft()` выпоняет прямое fft над последовательностями целых чисел от 1 до 8, 9, 25, 30, затем выполняется обратное преобразование над полученными результатами.
```text
длина 2^3
36, 10.4525, 5.65685, 4.32957, 4, 4.32957, 5.65685, 10.4525, // после прямого fft (выводится модуль комплексного числа)
1, 2, 3, 4, 5, 6, 7, 8,   //после обратного
длина 3^2
45, 13.1571, 7.00076, 5.19615, 4.56942, 4.56942, 5.19615, 7.00076, 13.1571, 
1, 2, 3, 4, 5, 6, 7, 8, 9, 
длина 5^2
325, 99.7341, 50.2634, 33.9559, 25.9469, 21.2663, 18.2602, 16.223, 14.8047, 13.8148, 13.1433, 12.7254, 12.5247, 12.5247, 12.7254, 13.1433, 13.8148, 14.8047, 16.223, 18.2602, 21.2663, 25.9469, 33.9559, 50.2634, 99.7341, 
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 
длина 2 * 3 * 5
465, 143.502, 72.146, 48.541, 36.8789, 30, 25.5195, 22.4171, 20.1845, 18.541, 17.3205, 16.4195, 15.7719, 15.3351, 15.0826, 15, 15.0826, 15.3351, 15.7719, 16.4195, 17.3205, 18.541, 20.1845, 22.4171, 25.5195, 30, 36.8789, 48.541, 72.146, 143.502
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
```
