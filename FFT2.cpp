#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
using std::cin, std::cout;

int reversebits(int x, int pw2) {
    int res = 0;
    for (int i = 0; i < pw2; i++) {
        if ((x >> i) & 1)
            res |= (1 << (pw2 - i - 1));
    }
    return res;
}

class Complex {
   public:
    long double r;
    long double i;

    Complex() : r(0), i(0) {}

    Complex(long double x) : r(x), i(0) {}

    Complex(long double r, long double i) : r(r), i(i) {}

    Complex(const Complex &tocopy) : r(tocopy.r), i(tocopy.i) {}

    Complex operator+(const Complex &second) const {
        return Complex(r + second.r, i + second.i);
    }

    Complex operator-(const Complex &second) const {
        return Complex(r - second.r, i - second.i);
    }

    Complex operator*(const Complex &second) const {
        return Complex(r * second.r - i * second.i, i * second.r + r * second.i);
    }

    Complex &operator*=(const Complex &second) {
        long double temp = r;
        long double temp1 = second.r;
        r = r * second.r - i * second.i;
        i = i * temp1 + temp * second.i;
        return *this;
    }

    Complex &operator=(const Complex &second) {
        r = second.r;
        i = second.i;
        return *this;
    }
};

std::ostream &operator<<(std::ostream &output, const Complex &val) {
    output << ((double)(int)(val.r * 100.0)) / 100 << "+" << ((double)(int)(val.i * 100.0)) / 100 << "i";
    return output;
}

void FFT(Complex *a, int n, Complex q) {
    int pw2 = 0;
    while ((1 << pw2) < n)
        pw2++;
    for (int i = 0; i < n; i++) {
        int rev = reversebits(i, pw2);
        if (i < rev)
            std::swap(a[i], a[rev]);
    }
    for (int l = 2; l <= n; l *= 2) {
        Complex cur = q;
        for (int ll = n; ll > l; ll /= 2)
            cur *= cur;
        for (int start = 0; start < n; start += l) {
            int mid = start + l / 2;
            Complex qdeg = 1;
            int pos = start;
            while (pos < mid) {
                Complex u = Complex(a[pos]);
                Complex v = Complex(a[pos + l / 2]) * qdeg;
                a[pos] = u + v;
                a[pos + l / 2] = u - v;
                qdeg *= cur;
                pos++;
            }
        }
    }
}

long long *mul(long long *p1, long long *p2, int sz) {
    int n = 1;
    while (n < sz)
        n *= 2;
    n *= 2;
    Complex *a = new Complex[n]();
    Complex *b = new Complex[n]();
    for (size_t i = 0; i < sz; i++)
        a[i] = p1[i];
    for (size_t i = 0; i < sz; i++)
        b[i] = p2[i];
    long double phi = 2 * acos(-1) / static_cast<long double>(n);
    Complex q = Complex(cos(phi), sin(phi));
    FFT(a, n, q);
    FFT(b, n, q);
    for (int i = 0; i < n; i++) {
        a[i] *= b[i];
    }
    FFT(a, n, Complex(cos(-phi), sin(-phi)));
    long long *res = new long long[n];
    for (int i = 0; i < n; i++) {
        res[i] = round(a[i].r / static_cast<long double>(n));
    }
    delete[] a;
    delete[] b;
    return res;
}
const int MAX = 100000;
int main() {
    std::ios::sync_with_stdio(0);
    cin.tie(0);
    int n, m, v;
    cin >> n;
    long long *a = new long long[MAX + 1]();
    long long *ra = new long long[MAX + 1]();
    for (int i = 0; i < n; i++) {
        cin >> v;
        ++a[v];
        ++ra[MAX - v];
    }
    long long *b = mul(a, ra, MAX + 1);
    cin >> m;
    for (int i = 0; i < m; i++) {
        cin >> v;
        if (v <= MAX)
            cout << b[MAX + v] << "\n";
        else
            cout << 0 << "\n";
    }
    return 0;
}