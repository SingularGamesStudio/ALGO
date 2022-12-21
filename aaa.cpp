#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
using std::cin, std::cout, std::complex;

int reversebits(int x, int pw2) {
    int res = 0;
    for (int i = 0; i < pw2; i++) {
        if ((x >> i) & 1)
            res |= (1 << (pw2 - i - 1));
    }
    return res;
}

long long mod = 7340033;

void FFT(complex<long double> *a, int n, complex<long double> q) {
    int pw2 = 0;
    while ((1 << pw2) < n)
        pw2++;
    for (int i = 0; i < n; i++) {
        int rev = reversebits(i, pw2);
        if (i < rev)
            std::swap(a[i], a[rev]);
    }
    for (int l = 2; l <= n; l *= 2) {
        complex<long double> cur = q;
        for (int ll = n; ll > l; ll /= 2)
            cur *= cur;
        for (int start = 0; start < n; start += l) {
            int mid = start + l / 2;
            complex<long double> qdeg = 1;
            int pos = start;
            while (pos < mid) {
                complex<long double> u = complex<long double>(a[pos]);
                complex<long double> v = complex<long double>(a[pos + l / 2]) * qdeg;
                a[pos] = u + v;
                a[pos + l / 2] = u - v;
                qdeg *= cur;
                pos++;
            }
        }
    }
}
const long long logN = 18;
long long *mul(long long *p1, long long *p2, int sz1, int sz2) {
    int n = (1ll << logN);
    complex<long double> *a = new complex<long double>[n]();
    complex<long double> *b = new complex<long double>[n]();
    for (size_t i = 0; i < sz1; i++)
        a[i] = p1[i];
    for (size_t i = 0; i < sz2; i++)
        b[i] = p2[i];
    long double phi = 2 * acos(-1) / static_cast<long double>(n);
    complex<long double> q = complex<long double>(cos(phi), sin(phi));
    FFT(a, n, q);
    FFT(b, n, q);
    for (int i = 0; i < n; i++) {
        a[i] *= b[i];
    }
    FFT(a, n, complex<long double>(cos(-phi), sin(-phi)));
    long long *res = new long long[n];
    for (int i = 0; i < n; i++) {
        res[i] = (static_cast<long long>(round(a[i].real() / static_cast<long double>(n))) % mod + mod) % mod;
    }
    delete[] a;
    delete[] b;
    return res;
}

long long gcd(long long a, long long b, long long &v1, long long &v2) {
    if (!a) {
        v1 = 0;
        v2 = 1;
        return b;
    }
    long long d = gcd(b % a, a, v2, v1);
    v1 -= (b / a) * v2;
    return d;
}

long long *reverse(long long *b, int n) {
    if (n == 1) {
        long long *rb = new long long[1];
        long long _;
        gcd(b[0], mod, rb[0], _);
        rb[0] = (rb[0] % mod + mod) % mod;
        return rb;
    }
    long long *b0 = new long long[n / 2];
    long long *b1 = new long long[n / 2];
    for (int i = 0; i < n / 2; i++) {
        b0[i] = b[i];
        b1[i] = b[n / 2 + i];
    }
    long long *b0r = reverse(b0, n / 2);
    long long *q = mul(b0, b0r, n / 2, n / 2);
    long long *b1b0r = mul(b1, b0r, n / 2, n / 2);
    for (int i = 0; i < n / 2; i++) {
        b1b0r[i] += q[n / 2 + i];
        b1b0r[i] %= mod;
    }
    long long *res = mul(b0r, b1b0r, n / 2, n / 2);
    for (int i = 0; i < n / 2; i++) {
        res[n / 2 + i] = (-res[i] + mod) % mod;
        res[i] = b0r[i];
    }
    delete[] b0;
    delete[] b1;
    delete[] b0r;
    delete[] b1b0r;
    delete[] q;
    return res;
}

int main() {
    std::ios::sync_with_stdio(0);
    cin.tie(0);
    int n, m, n0 = 1;  //
    cin >> m >> n;
    n++;
    while (n0 < n)
        n0 *= 2;
    while (n0 < m)
        n0 *= 2;
    n0 *= 2;
    long long *a = new long long[n0]();
    for (int i = 0; i < n; i++)
        cin >> a[i];
    if (a[0] == 0) {
        cout << "The ears of a dead donkey";
        return 0;
    }
    long long *ar = reverse(a, n0);
    for (int i = 0; i < m; i++) {
        cout << ar[i] << " ";
    }
    delete[] a;
    delete[] ar;
    return 0;
}