#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f(double x) {
    return log10(x + 1) - x + 0.5;
}

double df(double x) {
    return 1 / ((x + 1) * log(10)) - 1;
}

double d2f(double x) {
    return -1 / (log(10) * (x + 1) * (x + 1));
}

pair<double, int> newton_method(double a, double b, double eps) {
    double x0 = a;
    if (f(x0) * d2f(x0) <= 0) {
        x0 = b;
    }
    double x_prev = x0;
    double x = x_prev - f(x_prev) / df(x_prev);
    int k = 1;
    while (abs(x - x_prev) >= eps) {
        x_prev = x;
        x = x_prev - f(x_prev) / df(x_prev);
        ++k;
    }
    return {x, k};
}

double phi(double x) {
    return 0.5 + log10(x + 1);
}

// монотонна на [l, r]
double dphi(double x) {
    return 1 / ((x + 1) * log(10));
}

double golden_ratio(double a, double b, double eps, double (*f)(double)) {
    int k = 0;
    double y = a + (3 - sqrt(5)) / 2 * (b - a);
    double z = a + b - y;
    double f_y = (-1) * f(y);
    double f_z = (-1) * f(z);
    while (abs(a - b) > eps) {
        if (f_y <= f_z) {
            b = z;
            z = y;
            y = a + b - y;
            f_z = f_y;
            f_y = (-1) * f(y);
        } else {
            a = y;
            y = z;
            z = a + b - z;
            f_y = f_z;
            f_z = (-1) * f(z);
        }
        ++k;
    }
    return (a + b) / 2;
}

pair<double, int> simple_iterations(double l, double r, double eps) {
    double x_prev = (r - l) / 2;
    double x = phi(x_prev);
    int k = 1;
    // double q = max(dphi(l), dphi(r));
    double q = dphi(golden_ratio(l, r, eps, dphi));
    // q / (1 - q) * |x - x_prev| > eps
    while (q / (1 - q) * abs(x - x_prev) > eps) {
        x_prev = x;
        x = phi(x_prev);
        ++k;
    }
    return {x, k};
}

int main() {
    double eps;
    cin >> eps;
    auto [x1, k1] = newton_method(0.5, 0.9, eps);
    cout << "Solve newton_method: x = " << x1 << ", k = " << k1 << "\n";
    auto [x2, k2] = simple_iterations(0.5, 0.9, eps);
    cout << "Solve simple_iterations: x = " << x2 << ", k = " << k2 << "\n";
    return 0;
}