#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

double calc_exact_y(double x) {
    return 5 * x + x * x + x * x * log(abs(x));
}

double method_runge_romberg(double y1, double y2, int64_t p) {
    return (y1 - y2) / (pow(2, p) - 1);
}

vector<double> f(double x, double y, double z) {
    vector<double> F(2);
    F[0] = z;
    F[1] = (5 * x - 4 * y - 3 * x * z) / (x * x);
    return F;
}

void print(string method, vector<double> &X_h1, vector<double> &Y_h1, vector<double> &Y_h2, int64_t p) {
    cout << method << "\n"
         << "      x      |" << "      y\t   |" << "   exact y   |" << "\teps    | runge-romberg\n-------------------------------------------------------------------------\n";
    for (int i = 0; i < X_h1.size(); ++i) {
        double exact_y = calc_exact_y(X_h1[i]);
        cout << fixed << setprecision(9) << " " << X_h1[i] << " | " << Y_h1[i] << " | " << exact_y << " | " << abs(exact_y - Y_h1[i]) << " | " << method_runge_romberg(Y_h1[i], Y_h2[2 * i], p) << endl;
    }
}
