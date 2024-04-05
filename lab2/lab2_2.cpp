#include "../lab1/matrix.h"

using namespace std;

double f1(double x1, double x2) {
    return 2 * x1 * x1 - x2 + x2 * x2 - 2;
}

double df1_x1(double x1, double x2) {
    return 4 * x1;
}

double df1_x2(double x1, double x2) {
    return -1 + 2 * x2;
}

double f2(double x1, double x2) {
    return x1 - sqrt(x2 + 2) + 1;
}

double df2_x1(double x1, double x2) {
    return 1;
}

double df2_x2(double x1, double x2) {
    return -1 / (2 * sqrt(x2 + 2));
}

double norm(Matrix X) {
    double result = 0;
    for (int i = 0; i < X.GetRows(); ++i) {
        result += X(i, 0) * X(i, 0);
    }
    return sqrt(result);
}

pair<Matrix, int> newton_system(Matrix X, double eps) {
    int n = X.GetRows();
    Matrix J(n, n), B(n, 1), Delta_X(n, 1), X_Prev(n, 1);
    int k = 0;
    do {
        J(0, 0) = df1_x1(X(0, 0), X(1, 0));
        J(0, 1) = df1_x2(X(0, 0), X(1, 0));
        J(1, 0) = df2_x1(X(0, 0), X(1, 0));
        J(1, 1) = df2_x2(X(0, 0), X(1, 0));
        B(0, 0) = (-1) * f1(X(0, 0), X(1, 0));
        B(1, 0) = (-1) * f2(X(0, 0), X(1, 0));
        Delta_X = J.Solve(B);
        X_Prev = X;
        X = X + Delta_X;
        ++k;
    } while (norm(X - X_Prev) >= eps);
    return {X, k};
}

double phi1(double x1, double x2) {
    return sqrt(x2 + 2) - 1;
}

double phi2(double x1, double x2) {
    return 2 * x1 * x1 + x2 * x2 - 2;
}

double phi1_x1(double x1) {
    return 0;
}

double phi1_x2(double x2) {
    return 1 / (2 * sqrt(x2 + 2));
}

double phi2_x1(double x1) {
    return 4 * x1;
}

double phi2_x2(double x2) {
    return 2 * x2;
}

// double phi1(double x1, double x2) {
//     return sqrt(1 / 2 * (x2 - x2 * x2 + 2));
// }

// double phi2(double x1, double x2) {
//     return (x1 + 1) * (x1 + 1) - 2;
// }

// double phi1_x1(double x1) {
//     return 0;
// }

// double phi1_x2(double x2) {
//     return (1 - 2 * x2) / (4 * sqrt(1 / 2 * (x2 - x2 * x2 + 2)));
// }

// double phi2_x1(double x1) {
//     return 2 * (x1 + 1);
// }

// double phi2_x2(double x2) {
//     return 0;
// }

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

pair<Matrix, int> simple_iterations_system(double ax1, double bx1, double ax2, double bx2, Matrix X, double eps) {
    int n = X.GetRows();
    Matrix X_Prev(n, 1);
    int k = 0;
    Matrix M(n, n);
    M(0, 0) = phi1_x1(golden_ratio(ax1, bx1, eps, phi1_x1));
    M(0, 1) = phi1_x2(golden_ratio(ax2, bx2, eps, phi1_x2));
    M(1, 0) = phi2_x1(golden_ratio(ax1, bx1, eps, phi2_x1));
    M(1, 1) = phi2_x2(golden_ratio(ax2, bx2, eps, phi2_x2));
    // double q = M.norm();
    double q = 0.5;
    do {
        X_Prev = X;
        X(0, 0) = phi1(X_Prev(0, 0), X_Prev(1, 0));
        X(1, 0) = phi2(X_Prev(0, 0), X_Prev(1, 0));
        ++k;
    } while (q / (1 - q) * (X - X_Prev).norm() >= eps);
    return {X, k};
}

int main() {
    double eps;
    cin >> eps;
    Matrix X(2, 1);
    X(0, 0) = 1;
    X(1, 0) = 1;
    auto [X1, k1] = newton_system(X, eps);
    cout << "Solve newton_method:\n";
    cout << "X:\n"
         << X1 << "k = " << k1 << "\n";
    X(0, 0) = -0.3;
    X(1, 0) = -1.2;
    // auto [X2, k2] = simple_iterations_system(0.5, 1, 1, 1.5, X, eps);
    auto [X2, k2] = simple_iterations_system(-0.5, 0.5, -1, -0.5, X, eps);
    cout << "\nSolve simple_iterations:\n";
    cout << "X:\n"
         << X2 << "k = " << k2 << "\n";
    return 0;
}