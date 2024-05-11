#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

vector<double> f(double x, double y, double z) {
    vector<double> F(2);
    F[0] = z;
    F[1] = (5 * x - 4 * y - 3 * x * z) / (x * x);
    return F;
}

pair<vector<double>, vector<double>> method_euler(double x_start, double x_finish, double h, double y0, double z0) {
    int n = (x_finish - x_start) / h;
    vector<double> X(n + 1), Y(n + 1), Z(n + 1);
    for (int i = 0; i <= n; ++i) {
        X[i] = x_start + i * h;
    }
    Y[0] = y0;
    Z[0] = z0;
    for (int i = 1; i <= n; ++i) {
        vector<double> F = f(X[i - 1], Y[i - 1], Z[i - 1]);
        Y[i] = Y[i - 1] + h * F[0];
        Z[i] = Z[i - 1] + h * F[1];
    }
    return {X, Y};
}

int main() {
    auto [X, Y] = method_euler(1, 2, 0.1, 6, 8);
    for (int i = 0; i < X.size(); ++i) {
        cout << X[i] << "   " << Y[i] << "  " << 5 * X[i] + X[i] * X[i] + X[i] * X[i] * log(abs(X[i])) << endl;
    }
    return 0;
}