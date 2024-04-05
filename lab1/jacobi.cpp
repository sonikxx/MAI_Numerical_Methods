#include "./matrix.h"

using namespace std;

int main() {
    ifstream cin("test4.in");
    ofstream cout("test.out");
    double eps;
    cin >> eps;
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;
    auto [Res, k] = A.jacobi_method(eps);
    cout << "k = " << k << "\n";
    cout << "A:\n";
    cout << Res.first;
    // check Lambda = U_T * A * U
    cout << "\nU_T * A * U:\n";
    cout << Res.second.Transpose().MulMatrixReturn(A).MulMatrixReturn(Res.second);
    cout << "\nSelf_value:\n";
    for (int i = 0; i < Res.first.GetRows(); ++i) {
        cout << Res.first(i, i) << ' ';
    }
    cout << "\n";
    cout << "\nSelf_vectors:\n";
    cout << Res.second;
    return 0;
}