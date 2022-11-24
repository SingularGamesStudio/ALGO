#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using std::cin, std::cout;

const long long MOD = 998244353;
const long long logN = 21;
const long long root1 = 421;
const long long revroot1 = 884430270;

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

    long long divN(long long n) {
        return static_cast<long long>(r / static_cast<long double>(n) + static_cast<long double>(0.5));
    }

    static Complex get1root() {
        long double phi = 2 * acos(-1) / static_cast<long double>(1ll << logN);
        return Complex(cos(phi), sin(phi));
    }

    static Complex get1rootrev() {
        long double phi = -2 * acos(-1) / static_cast<long double>(1ll << logN);
        return Complex(cos(phi), sin(phi));
    }
};

std::ostream &operator<<(std::ostream &output, const Complex &val) {
    output << ((double)(int)(val.r * 100.0)) / 100 << "+" << ((double)(int)(val.i * 100.0)) / 100 << "i";
    return output;
}

long long gcd(long long a, long long b, long long &x, long long &y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    long long x1, y1;
    long long d = gcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}

class Zmod {
   public:
    long long val;

    Zmod() : val(0) {}

    Zmod(long long x) : val(x) {}

    Zmod(const Zmod &tocopy) : val(tocopy.val) {}

    Zmod operator+(const Zmod &second) const {
        return Zmod((val + second.val) % MOD);
    }

    Zmod operator-(const Zmod &second) const {
        return Zmod((val - second.val + MOD) % MOD);
    }

    Zmod operator*(const Zmod &second) const {
        return Zmod((val * second.val) % MOD);
    }

    Zmod &operator*=(const Zmod &second) {
        val = (val * second.val) % MOD;
        return *this;
    }

    Zmod &operator=(const Zmod &second) {
        val = second.val;
        return *this;
    }

    long long divN(long long n) {
        long long x, y;
        gcd(n, MOD, x, y);
        long long revn = (x % MOD + MOD) % MOD;
        return (val * revn) % MOD;
    }

    static Zmod get1root() {
        return Zmod(root1);
    }

    static Zmod get1rootrev() {
        return Zmod(revroot1);
    }
};
template <typename T>
void FFT(T *a, int n, T q) {
    int pw2 = 0;
    while ((1 << pw2) < n)
        pw2++;
    for (int i = 0; i < n; i++) {
        int rev = reversebits(i, pw2);
        if (i < rev)
            std::swap(a[i], a[rev]);
    }
    for (int l = 2; l <= n; l *= 2) {
        T cur = q;
        for (int ll = n; ll > l; ll /= 2)
            cur *= cur;
        for (int start = 0; start < n; start += l) {
            int mid = start + l / 2;
            T qdeg = 1;
            int pos = start;
            while (pos < mid) {
                T u = T(a[pos]);
                T v = T(a[pos + l / 2]) * qdeg;
                a[pos] = u + v;
                a[pos + l / 2] = u - v;
                qdeg *= cur;
                pos++;
            }
        }
    }
}

template <typename T>
long long *mul(long long *p1, long long *p2, int sz) {
    int n = (1ll << logN);
    T *a = new T[n]();
    T *b = new T[n]();
    for (size_t i = 0; i < sz; i++)
        a[i] = p1[i];
    for (size_t i = 0; i < sz; i++)
        b[i] = p2[i];
    T q = T::get1root();
    FFT(a, n, q);
    FFT(b, n, q);
    for (int i = 0; i < n; i++) {
        a[i] *= b[i];
    }
    q = T::get1rootrev();
    FFT(a, n, q);
    long long *res = new long long[n];
    for (int i = 0; i < n; i++) {
        res[i] = a[i].divN(n);
    }
    delete[] a;
    delete[] b;
    return res;
}

int binpow(long long a, int p, long long mod) {
    if (p == 0)
        return 1 % mod;
    long long k = binpow(a, p / 2, mod);
    k = (k * k) % mod;
    if (p % 2)
        k = (k * a) % mod;
    return k;
}

long long stress(int n, int m) {
    long long ans = 0;
    for (int x = 1; x < m; x++) {
        for (int y = x; y < m; y++) {
            for (int z = 1; z < m; z++) {
                if ((binpow(x, n, m) + binpow(y, n, m)) % m == binpow(z, n, m))
                    ans++;
            }
        }
    }
    return ans;
}

/*int main() {
    for (long long i = 2; i < MOD; i++) {
        if ((i * root1) % MOD == 1ll) {
            cout << i;
            return 0;
        }
    }
}*/

int main() {
    std::ios::sync_with_stdio(0);
    cin.tie(0);
    std::mt19937 rnd(42);
    // while (1) {
    int n = rnd() % 10, m = rnd() % 1000;
    cin >> n >> m;
    long long *a = new long long[m]();
    for (int i = 1; i < m; i++) {
        ++a[binpow(i, n, m)];
    }
    long long *b = mul<Zmod>(a, a, m);
    long long ans = 0, ans1 = 0;
    for (int i = 0; i < m; i++) {
        // cout << a[i] << " (" << b[i] << ", " << b[i + m] << ")\n";
        ans += a[i] * (b[i] + b[i + m]);
        ans1 += a[((i * 2) % m)] * a[i];
    }
    long long res = (ans - ans1) / 2ll + ans1;
    cout << res;
    /*if (res != stress(n, m)) {
        cout << "stupid: " << stress(n, m) << "\nclever: " << res << "\n"
             << n << " " << m;
        return 0;
    }*/
    //}
    return 0;
}
