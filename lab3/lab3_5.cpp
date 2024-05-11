#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f(double x) {
    return (x * x * sqrt(36 - x * x));
}

double rectangle_method(double a, double b, double step) {
    double res = 0;
    double x_prev = a;
    // F = sum(h * f((x_i-1 + x_i)/2))
    for (double cur = a + step; cur <= b; cur += step) {
        res += step * f((x_prev + cur) / 2);
        x_prev = cur;
    }
    return res;
}

double trapeze_method(double a, double b, double step) {
    double res = 0;
    double x_prev = a;
    // F = 1/2 * sum(h * (f(x_i) + f(x_i-1)))
    for (double cur = a + step; cur <= b; cur += step) {
        res += step * (f(cur) + f(x_prev));
        x_prev = cur;
    }
    return 0.5 * res;
}

double simpson_method(double a, double b, double step) {
    double res = 0;
    // F = 1/3 * sum(h * (f(x_i-1) + 4 * f(x_i-1/2) + f(x_i)))
    for (double cur = a + 2 * step; cur <= b; cur += 2 * step) {
        res += step * (f(cur - 2 * step) + 4 * f(cur - step) + f(cur));
    }
    return (1.0 / 3 * res);
}

// p - порядок погрешности
double runge_romberg_method(double F1, double F2, double h1, double h2, double p) {
    // F = F1 + (F1 - F2) / (k^p - 1)
    return (F1 + (F1 - F2) / (pow(h2 / h1, p) - 1));
}

int main() {
    double a = 1;
    double b = 5;
    double h1 = 1.0, h2 = 0.5;
    double absolute_value = 186.625164060908;
    cout << "Rectangle method:\n";
    double F1 = rectangle_method(a, b, h1);
    double F2 = rectangle_method(a, b, h2);
    double F = runge_romberg_method(F1, F2, h1, h2, 2);
    cout << "F = " << F1 << " with h = " << h1 << "\n";
    cout << "F = " << F2 << " with h = " << h2 << "\n";
    cout << "error = " << abs(F - absolute_value) << "\n";
    cout << "\nTrapeze method:\n";
    F1 = trapeze_method(a, b, h1);
    F2 = trapeze_method(a, b, h2);
    F = runge_romberg_method(F1, F2, h1, h2, 2);
    cout << "F = " << F1 << " with h = " << h1 << "\n";
    cout << "F = " << F2 << " with h = " << h2 << "\n";
    cout << "error = " << abs(F - absolute_value) << "\n";
    cout << "\nSimpson method:\n";
    F1 = simpson_method(a, b, h1);
    F2 = simpson_method(a, b, h2);
    F = runge_romberg_method(F1, F2, h1, h2, 4);
    cout << "F = " << F1 << " with h = " << h1 << "\n";
    cout << "F = " << F2 << " with h = " << h2 << "\n";
    cout << "error = " << abs(F - absolute_value) << "\n";
    return 0;
}