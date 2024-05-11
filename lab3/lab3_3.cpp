#include "../../../matplotlibcpp.h"
#include "../lab1/matrix.h"
#include <cmath>
#include <iostream>
#include <vector>

namespace plt = matplotlibcpp;
using namespace std;

vector<double> mnk(vector<double> &x, vector<double> &y, int m) {
    int n = x.size();
    ++m; // для построения многочлена степени m - нужна m+1 точка
    Matrix Y(n, 1), Phi(n, m);
    for (int i = 0; i < n; ++i) {
        Y(i, 0) = y[i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            Phi(i, j) = pow(x[i], j);
        }
    }
    Matrix Phi_T = Phi.Transpose();
    // (Ф_T * Ф) * a = Ф_T * Y
    Matrix B = Phi_T.MulMatrixReturn(Y);
    Matrix res = (Phi_T.MulMatrixReturn(Phi)).Solve(B);
    vector<double> polynom;
    for (int i = 0; i < res.GetRows(); ++i) {
        for (int j = 0; j < res.GetCols(); ++j) {
            polynom.push_back(res(i, j));
        }
    }
    return polynom;
}

vector<double> at(vector<double> &x, vector<double> &p) {
    vector<double> res;
    for (int i = 0; i < x.size(); ++i) {
        double tmp = 0;
        for (int j = 0; j < p.size(); ++j) {
            tmp += p[j] * pow(x[i], j);
        }
        res.push_back(tmp);
    }
    return res;
}

void error(vector<double> &x, vector<double> &y, vector<double> &p) {
    double f_x = 0, error = 0;
    for (int i = 0; i < x.size(); ++i) {
        f_x = 0;
        for (int j = 0; j < p.size(); ++j) {
            f_x += p[j] * pow(x[i], j);
        }
        error += (f_x - y[i]) * (f_x - y[i]);
    }
    cout << "error = " << error << "\n\n";
}

int main() {
    vector<double> x = {0.0, 1.7, 3.4, 5.1, 6.8, 8.5};
    vector<double> split_x;
    for (double i = x[0]; i <= x[x.size() - 1]; i += 0.1) {
        split_x.push_back(i);
    }
    vector<double> y = {0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155};
    vector<double> polynom = mnk(x, y, 1);
    cout << "MNK polynom 1:\n";
    for (int i = 0; i < polynom.size(); ++i) {
        if (i != polynom.size() - 1)
            cout << polynom[i] << "x^" << i << " + ";
        else
            cout << polynom[i] << "x^" << i << "\n";
    }
    error(x, y, polynom);
    vector<double> z = at(split_x, polynom);
    polynom.clear();
    polynom = mnk(x, y, 2);
    cout << "MNK polynom 2:\n";
    for (int i = 0; i < polynom.size(); ++i) {
        if (i != polynom.size() - 1)
            cout << polynom[i] << "x^" << i << " + ";
        else
            cout << polynom[i] << "x^" << i << "\n";
    }
    error(x, y, polynom);
    vector<double> w = at(split_x, polynom);
    plt::scatter(x, y, 10, {{"label", "f(x)"}, {"color", "red"}});
    plt::plot(split_x, z, {{"label", "polynom 1"}});
    plt::plot(split_x, w, {{"label", "polynom 2"}});
    plt::legend();
    plt::show();
    return 0;
}

// g++ lab3_3.cpp -I//usr/include/python3.10 -lpython3.10 && ./a.out
