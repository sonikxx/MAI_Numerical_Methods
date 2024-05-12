#include "./lab4_1.h"

vector<double> num_vector(vector<double> a, double n) {
    for (int i = 0; i < a.size(); ++i) {
        a[i] *= n;
    }
    return a;
}

vector<vector<double>> method_runge_kutta_4(double x_start, double x_finish, double h, double y0, double z0, int iter = -1) {
    int n;
    if (iter == -1) {
        n = (x_finish - x_start) / h;
    } else
        n = iter;
    vector<double> X(n + 1), Y(n + 1), Z(n + 1);
    for (int i = 0; i <= n; ++i) {
        X[i] = x_start + i * h;
    }
    Y[0] = y0;
    Z[0] = z0;
    for (int i = 1; i <= n; ++i) {
        vector<double> K_1 = num_vector(f(X[i - 1], Y[i - 1], Z[i - 1]), h);
        vector<double> K_2 = num_vector(f(X[i - 1] + h / 2, Y[i - 1] + K_1[0] / 2, Z[i - 1] + K_1[1] / 2), h);
        vector<double> K_3 = num_vector(f(X[i - 1] + h / 2, Y[i - 1] + K_2[0] / 2, Z[i - 1] + K_2[1] / 2), h);
        vector<double> K_4 = num_vector(f(X[i - 1] + h, Y[i - 1] + K_3[0], Z[i - 1] + K_3[1]), h);
        Y[i] = Y[i - 1] + 1.0 / 6 * (K_1[0] + 2 * K_2[0] + 2 * K_3[0] + K_4[0]);
        Z[i] = Z[i - 1] + 1.0 / 6 * (K_1[1] + 2 * K_2[1] + 2 * K_3[1] + K_4[1]);
    }
    return {X, Y, Z};
}

pair<vector<double>, vector<double>> method_adams(double x_start, double x_finish, double h, double y0, double z0) {
    int n = (x_finish - x_start) / h;
    vector<double> X(n + 1), Y(n + 1), Z(n + 1);
    for (int i = 0; i <= n; ++i) {
        X[i] = x_start + i * h;
    }
    vector<vector<double>> res_runge_kutta = method_runge_kutta_4(x_start, x_finish, h, y0, z0, 4);
    Y[0] = res_runge_kutta[1][0];
    Y[1] = res_runge_kutta[1][1];
    Y[2] = res_runge_kutta[1][2];
    Y[3] = res_runge_kutta[1][3];
    Z[0] = res_runge_kutta[2][0];
    Z[1] = res_runge_kutta[2][1];
    Z[2] = res_runge_kutta[2][2];
    Z[3] = res_runge_kutta[2][3];
    for (int i = 4; i <= n; ++i) {
        vector<double> F_k = f(X[i - 1], Y[i - 1], Z[i - 1]);
        vector<double> F_k_1 = f(X[i - 2], Y[i - 2], Z[i - 2]);
        vector<double> F_k_2 = f(X[i - 3], Y[i - 3], Z[i - 3]);
        vector<double> F_k_3 = f(X[i - 4], Y[i - 4], Z[i - 4]);
        Y[i] = Y[i - 1] + h / 24 * (55 * F_k[0] - 59 * F_k_1[0] + 37 * F_k_2[0] - 9 * F_k_3[0]);
        Z[i] = Z[i - 1] + h / 24 * (55 * F_k[1] - 59 * F_k_1[1] + 37 * F_k_2[1] - 9 * F_k_3[1]);
    }
    return {X, Y};
}

int main() {
    vector<vector<double>> res_h1 = method_runge_kutta_4(1, 2, 0.1, 6, 8);
    vector<double> X = res_h1[0];
    vector<double> Y_h1 = res_h1[1];
    vector<vector<double>> res_h2 = method_runge_kutta_4(1, 2, 0.05, 6, 8);
    vector<double> Y_h2 = res_h2[1];
    print("Method Runge-Kutta:", X, Y_h1, Y_h2, 4);
    auto [X2_h1, Y2_h1] = method_adams(1, 2, 0.1, 6, 8);
    auto [X2_h2, Y2_h2] = method_adams(1, 2, 0.05, 6, 8);
    print("\nMethod Adams:", X2_h1, Y2_h1, Y2_h2, 4);
    return 0;
}