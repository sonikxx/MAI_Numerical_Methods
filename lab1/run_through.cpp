#include "./matrix.h"

using namespace std;

int main() {
    ifstream cin("test2.in");
    ofstream cout("test.out");
    int n;
    cin >> n;
    Matrix A(n, n);
    cin >> A;
    Matrix B(n, 1);
    cin >> B;
    cout << A.run_through_method(B);
    return 0;
}