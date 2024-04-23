#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double df(vector<double> &x, vector<double> &y, double dot) {
    auto it = lower_bound(x.begin(), x.end(), dot);
    int interval = it - x.begin() - 1;
    if (interval == -1)
        interval = 0;
    double a = (y[interval + 1] - y[interval]) / (x[interval + 1] - x[interval]);
    double b = ((y[interval + 2] - y[interval + 1]) / (x[interval + 2] - x[interval + 1]) - (y[interval + 1] - y[interval]) / (x[interval + 1] - x[interval])) / (x[interval + 2] - x[interval]);
    return (a + b * (2 * dot - x[interval] - x[interval + 1]));
}

double d2f(vector<double> &x, vector<double> &y, double dot) {
    auto it = lower_bound(x.begin(), x.end(), dot);
    int interval = it - x.begin() - 1;
    if (interval == -1)
        interval = 0;
    // y``(x) = 2 * ((y_i+2 - y_i+1)/(x_i+2 - x_i+1) - (y_i+1 - y_i)/(x_i+1 - x_i)) / (x_i+2 - x_i)
    return (2 * (((y[interval + 2] - y[interval + 1]) / (x[interval + 2] - x[interval + 1])) - ((y[interval + 1] - y[interval])) / (x[interval + 1] - x[interval]))) / (x[interval + 2] - x[interval]);
}

int main() {
    vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    vector<double> y = {0.0, 0.86603, 1.0, 0.0, -2.0};
    double dot = 2.0;
    cout << "y`(x) = " << df(x, y, dot) << "\n";
    cout << "y``(x) = " << d2f(x, y, dot) << "\n";
    return 0;
}