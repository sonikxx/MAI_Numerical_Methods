#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

vector<double> mult(vector<double> a, vector<double> b) {
    vector<double> res(a.size() + b.size() - 1);
    for (int i = 0; i < b.size(); ++i) {
        for (int j = 0; j < a.size(); ++j) {
            res[i + j] += a[j] * b[i];
        }
    }
    return res;
}

vector<double> mult_number(vector<double> a, double n) {
    vector<double> res(a.size());
    for (int i = 0; i < a.size(); ++i) {
        res[i] = a[i] * n;
    }
    return res;
}

vector<double> sum(vector<double> a, vector<double> b) {
    int n_max = max(a.size(), b.size());
    int n_min = min(a.size(), b.size());
    vector<double> res(n_max);
    for (int i = 0; i < n_min; ++i) {
        res[i] = a[i] + b[i];
    }
    for (int i = n_min; i < n_max; ++i) {
        if (a.size() > b.size()) {
            res[i] = a[i];
        } else {
            res[i] = b[i];
        }
    }
    return res;
}

vector<double> lagrange_polynom(vector<double> &x, vector<double> &y) {
    int n = x.size();
    vector<double> L(n);
    for (int i = 0; i < n; ++i) {
        vector<double> li;
        li.push_back(1);
        // l0 = (x - x1) / (x0 - x1) * (x - x2) / (x1 - x2) * ...
        double c = 1; // общий знаменатель
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                li = mult(li, {(-1) * x[j], 1});
                c *= (x[i] - x[j]);
            }
        }
        L = sum(L, mult_number(li, y[i] / c));
    }
    return L;
}

vector<double> newton_polynom(vector<double> &x, vector<double> &y) {
    int n = x.size();
    vector<vector<double>> table(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        table[i][0] = y[i];
    }
    // f(xi, xj) = (fi - fj) / (xi - xj)
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / (x[i] - x[i + j]);
        }
    }
    vector<double> P(n);
    vector<double> k;
    // P = f(x0) + (x - x0) * f(x0, x1) + (x - x0) * (x - x1) * f(x0, x1, x2) + ...
    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            k.push_back(1);
        } else {
            k = mult(k, {(-1) * x[i - 1], 1});
        }
        P = sum(P, mult_number(k, table[0][i]));
    }
    return P;
}

double f(double x) {
    return 1 / (x * x) + x * x;
}

int main() {
    vector<double> x = {0.1, 0.5, 0.9, 1.3};
    vector<double> y;
    for (int i = 0; i < x.size(); ++i) {
        y.push_back(f(x[i]));
    }
    cout << "Polynom Lagrange:\n";
    vector<double> polynom_l = lagrange_polynom(x, y);
    for (int i = 0; i < polynom_l.size(); ++i) {
        if (i != polynom_l.size() - 1)
            cout << polynom_l[i] << "x^" << i << " + ";
        else
            cout << polynom_l[i] << "x^" << i << "\n";
    }
    double res = 0;
    double X = 0.8;
    for (int i = 0; i < polynom_l.size(); ++i) {
        res += polynom_l[i] * pow(X, i);
    }
    cout << "f(x*) = " << res << "\n";
    cout << "error = " << abs(f(X) - res) << "\n";

    cout << "\nPolynom Newton:\n";
    vector<double> polynom_n = newton_polynom(x, y);
    for (int i = 0; i < polynom_n.size(); ++i) {
        if (i != polynom_n.size() - 1)
            cout << polynom_n[i] << "x^" << i << " + ";
        else
            cout << polynom_n[i] << "x^" << i << "\n";
    }
    res = 0;
    for (int i = 0; i < polynom_n.size(); ++i) {
        res += polynom_n[i] * pow(X, i);
    }
    cout << "f(x*) = " << res << "\n";
    cout << "error = " << abs(f(X) - res) << "\n";
    return 0;
}