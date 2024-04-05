#include "./matrix.h"

using namespace std;

int main() {
    ifstream cin("test1.in");
    ofstream cout("test.out");
    string flag;
    cin >> flag;
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;
    if (flag == "solve") {
        Matrix B(n, 1);
        cin >> B;
        cout << A.Solve(B);
        auto [L, U] = A.LU();
        cout << "L:\n";
        cout << L;
        cout << "U:\n";
        cout << U;
        cout << "L*U:\n";
        L.MulMatrix(U);
        vector<int> swp = A.GetSwp();
        for (int i = 0; i < swp.size(); ++i) {
            swap(L(i), L(swp[i]));
        }
        cout << L;
        if (L == A) {
            cout << "success\n";
        } else {
            cout << "fail\n";
        }
    } else if (flag == "det") {
        cout << A.Determinant();
    } else if (flag == "inverse") {
        cout << A.InverseMatrix();
    }
    return 0;
}