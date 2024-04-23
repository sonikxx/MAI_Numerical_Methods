#include "../lab1/matrix.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double cubic_spline(vector<double> &x, vector<double> &y, double dot) {
    int n = x.size() - 1;
    vector<double> a(n), b(n), c(n), d(n), h(n);
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i]; // нумерация с 1
    }
    Matrix A(n - 1, n - 1), B(n - 1, 1);
    // A(i, i) = 2 * (hi + hi-1), A(i,i-1) = hi, A(i, i+1) = hi+1
    // Bi = 3 * ((fi - fi-1) / hi - (fi-1 - fi-2) / hi-1)
    for (int i = 0; i < n - 1; ++i) {
        A(i, i) = 2 * (h[i] + h[i + 1]);
        if (i > 0)
            A(i, i - 1) = h[i];
        if (i < n - 2)
            A(i, i + 1) = h[i + 1];
        B(i, 0) = 3 * ((y[i + 2] - y[i + 1]) / h[i + 1] - (y[i + 1] - y[i]) / h[i]);
    }
    Matrix C = A.run_through_method(B);
    for (int i = 0; i < n; ++i) {
        if (i == 0)
            c[i] = 0;
        else
            c[i] = C(i - 1, 0);
    }
    for (int i = 0; i < n; ++i) {
        a[i] = y[i];
        if (i < n - 1) {
            b[i] = (y[i + 1] - y[i]) / h[i] - 1.0 / 3 * h[i] * (c[i + 1] + 2 * c[i]);
            d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
        } else {
            b[i] = (y[i + 1] - y[i]) / h[i] - 2.0 / 3 * h[i] * c[i];
            d[i] = (-1) * c[i] / (3 * h[i]);
        }
    }
    auto it = lower_bound(x.begin(), x.end(), dot);
    int interval = it - x.begin() - 1;
    if (interval == -1)
        interval = 0;
    cout << "S(x) = " << a[interval] << " + " << b[interval] << " * (x - " << x[interval] << ") + " << c[interval] << " * (x - " << x[interval] << ")^2 + " << d[interval] << " * (x - " << x[interval] << ")^3\n";
    // S(x) = ai + bi * (x - xi) + ci * (x - xi)^2 + di * (x - xi)^3
    double res = a[interval] + b[interval] * (dot - x[interval]) + c[interval] * pow((dot - x[interval]), 2) + d[interval] * pow((dot - x[interval]), 3);
    return res;
}

int main() {
    vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7};
    vector<double> y = {100.01, 4.2500, 2.0446, 2.2817, 3.2360};
    double dot = 0.8;
    double res = cubic_spline(x, y, dot);
    cout << "y(x*) = " << res << "\n";
    return 0;
}