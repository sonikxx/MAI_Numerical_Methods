#include "./matrix.h"

using namespace std;

int main() {
    ifstream cin("test5.in");
    ofstream cout("test.out");
    double eps;
    cin >> eps;
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;
    auto [Q, R] = A.qr_decomposition();
    cout << "Q:\n";
    cout << Q;
    cout << "\nR:\n";
    cout << R;
    cout << "\nQ * R:\n";
    cout << Q.MulMatrixReturn(R);
    vector<complex<double>> labmda = A.qr_method(eps);
    cout << "\nLambda:\n";
    for (int i = 0; i < labmda.size(); ++i) {
        cout << labmda[i] << " ";
    }
    cout << "\n";
    return 0;
}