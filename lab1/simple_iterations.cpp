#include "./matrix.h"

using namespace std;

int main() {
    ifstream cin("test3.in");
    ofstream cout("test.out");
    double eps;
    cin >> eps;
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;
    Matrix B(n, 1);
    cin >> B;
    cout << "Simple iterations methods:\n";
    auto [Res1, k1] = A.simple_iterations(B, eps);
    cout << "k = " << k1 << "\n";
    cout << Res1;
    cout << "\nSeidel methods:\n";
    auto [Res2, k2] = A.seidel(B, eps);
    cout << "k = " << k2 << "\n";
    cout << Res2;
    return 0;
}