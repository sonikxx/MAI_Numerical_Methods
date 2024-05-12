#include "./lab4_1.h"

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

pair<vector<double>, vector<double>> first_improved_method_euler(double x_start, double x_finish, double h, double y0, double z0) {
    h = h / 2;
    int n = (x_finish - x_start) / h;
    vector<double> X(n + 1), Y(n + 1), Z(n + 1);
    for (int i = 0; i <= n; ++i) {
        X[i] = x_start + i * h;
    }
    Y[0] = y0;
    Z[0] = z0;
    for (int i = 1; i <= n; ++i) {
        vector<double> F = f(X[i - 1], Y[i - 1], Z[i - 1]);
        if (i % 2 == 0) {
            Y[i] = Y[i - 2] + h * 2 * F[0];
            Z[i] = Z[i - 2] + h * 2 * F[1];
        } else {
            Y[i] = Y[i - 1] + h * F[0];
            Z[i] = Z[i - 1] + h * F[1];
        }
    }
    return {X, Y};
}

pair<vector<double>, vector<double>> first_method_euler_cauchy(double x_start, double x_finish, double h, double y0, double z0) {
    int n = (x_finish - x_start) / h;
    vector<double> X(n + 1), Y(n + 1), Y_tmp(n + 1), Z(n + 1), Z_tmp(n + 1);
    for (int i = 0; i <= n; ++i) {
        X[i] = x_start + i * h;
    }
    Y[0] = y0;
    Z[0] = z0;
    for (int i = 1; i <= n; ++i) {
        vector<double> F_tmp = f(X[i - 1], Y[i - 1], Z[i - 1]);
        Y_tmp[i] = Y[i - 1] + h * F_tmp[0];
        Z_tmp[i] = Z[i - 1] + h * F_tmp[1];
        vector<double> F = f(X[i], Y_tmp[i], Z_tmp[i]);
        Y[i] = Y[i - 1] + h * (F_tmp[0] + F[0]) / 2;
        Z[i] = Z[i - 1] + h * (F_tmp[1] + F[1]) / 2;
    }
    return {X, Y};
}

int main() {
    auto [X_h1, Y_h1] = method_euler(1, 2, 0.1, 6, 8);
    auto [X_h2, Y_h2] = method_euler(1, 2, 0.05, 6, 8);
    print("Method Euler:", X_h1, Y_h1, Y_h2, 1);
    auto [X2_h1, Y2_h1] = first_improved_method_euler(1, 2, 0.1, 6, 8);
    auto [X2_h2, Y2_h2] = first_improved_method_euler(1, 2, 0.05, 6, 8);
    print("\n\nFirst improved method Euler:", X2_h1, Y2_h1, Y2_h2, 2);
    auto [X3_h1, Y3_h1] = first_method_euler_cauchy(1, 2, 0.1, 6, 8);
    auto [X3_h2, Y3_h2] = first_method_euler_cauchy(1, 2, 0.05, 6, 8);
    print("\n\nFirst method Euler-Cauchy:", X3_h1, Y3_h1, Y2_h2, 2);
    return 0;
}