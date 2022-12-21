#include <limits.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <compare>
#include <complex>
#include <deque>
#include <iostream>
#include <string>
#include <vector>

using std::vector, std::string, std::complex, std::cin, std::cout, std::deque;

enum signs { neg = -1, zero = 0, pos = 1 };

signs mulsigns(signs a, signs b) {
    if (a == signs::zero || b == signs::zero)
        return signs::zero;
    if (a == b)
        return signs::pos;
    return signs::neg;
}

int reversebits(int x, int pw2) {
    int res = 0;
    for (int i = 0; i < pw2; i++) {
        if ((x >> i) & 1)
            res |= (1 << (pw2 - i - 1));
    }
    return res;
}

int MOD = 128;
const int DIGITS = 7;

void FFT(complex<double> *a, int n, complex<double> q) {
    int pw2 = 0;
    while ((1 << pw2) < n)
        pw2++;
    for (int i = 0; i < n; i++) {
        int rev = reversebits(i, pw2);
        if (i < rev)
            std::swap(a[i], a[rev]);
    }
    for (int l = 2; l <= n; l *= 2) {
        complex<double> cur = q;
        for (int ll = n; ll > l; ll /= 2)
            cur *= cur;
        for (int start = 0; start < n; start += l) {
            int mid = start + l / 2;
            complex<double> qdeg = 1;
            int pos = start;
            while (pos < mid) {
                complex<double> u = complex<double>(a[pos]);
                complex<double> v = complex<double>(a[pos + l / 2]) * qdeg;
                a[pos] = u + v;
                a[pos + l / 2] = u - v;
                qdeg *= cur;
                pos++;
            }
        }
    }
}

class PoweredInteger;

class BigInteger {
private:
    int get(size_t id) const {
        if (id < data.size())
            return data[id];
        return 0;
    }

    void upgrade(int x, int &zeros) {
        if (x) {
            for (int j = 0; j < zeros; j++) {
                data.push_back(0);
            }
            data.push_back(x);
            zeros = 0;
        } else
            ++zeros;
    }

    void add1() {  // add 1 to the abs
        if (sign == signs::zero) {
            data.push_back(1);
            sign = signs::pos;
            return;
        }
        int delta = 1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            delta = data[i] / MOD;
            data[i] = data[i] % MOD;
            if (delta == 0)
                break;
        }
        if (delta != 0)
            data.push_back(delta);
    }

    void sub1() {  // substract 1 from the abs
        if (sign == signs::zero) {
            data.push_back(1);
            sign = signs::neg;
            return;
        }
        int delta = -1;
        for (size_t i = 0; i < data.size(); i++) {
            data[i] += delta;
            if (data[i] < 0) {
                delta = -1;
                data[i] += MOD;
            } else
                break;
        }
        while (data.size() && data[data.size() - 1] == 0)
            data.pop_back();
        data.shrink_to_fit();
        if (data.size() == 0)
            sign = signs::zero;
    }

    void swap(BigInteger &other) {
        std::swap(data, other.data);
        std::swap(sign, other.sign);
    }

    std::strong_ordering cmpabs(const BigInteger &second) const {
        if (data.size() != second.data.size())
            return data.size() <=> second.data.size();
        for (int i = data.size() - 1; i >= 0; --i) {
            if (data[i] != second.data[i])
                return data[i] <=> second.data[i];
        }
        return std::strong_ordering::equal;
    }

    BigInteger add(const BigInteger &second, signs sign1 = signs::pos) const {
        if (second.sign == signs::zero || sign1 == signs::zero)
            return BigInteger(*this);
        if (sign == signs::zero)
            return BigInteger(second);
        BigInteger res = BigInteger();
        int size = std::max(data.size(), second.data.size()) + 1;
        int delta = 0;
        int zeros = 0;
        signs sign2 = mulsigns(second.sign, sign1);
        if (sign == sign2) {  // sum of modules
            for (int i = 0; i < size; i++) {
                int now = get(i) + second.get(i) + delta;
                delta = now / MOD;
                now = now % MOD;
                res.upgrade(now, zeros);
            }
            res.sign = sign;
        } else {  // diff of modules
            const BigInteger *abs1 = this;
            const BigInteger *abs2 = &second;
            bool swapped = false;
            if (cmpabs(second) == std::strong_ordering::less) {
                std::swap(abs1, abs2);
                swapped = true;
            }
            for (int i = 0; i < size; i++) {
                int now = (*abs1).get(i) - (*abs2).get(i) + delta;
                if (now < 0) {
                    delta = -1;
                    now = now + MOD;
                } else
                    delta = 0;
                res.upgrade(now, zeros);
            }
            if (zeros == size) {
                res.sign = zero;
            } else
                res.sign = swapped ? sign2 : sign;
        }
        return res;
    }

    BigInteger mul(const BigInteger &second) const {
        if (sign == signs::zero || second.sign == signs::zero)
            return 0;
        size_t n = 1;
        while (n < data.size() || n < second.data.size())
            n *= 2;
        n *= 2;
        complex<double> *a = new complex<double>[n]();
        complex<double> *b = new complex<double>[n]();
        for (size_t i = 0; i < data.size(); i++)
            a[i] = data[i];
        for (size_t i = 0; i < second.data.size(); i++)
            b[i] = second.data[i];
        double phi = 2 * acos(-1) / static_cast<double>(n);
        complex<double> q = complex<double>(cos(phi), sin(phi));
        FFT(a, n, q);
        FFT(b, n, q);
        for (size_t i = 0; i < n; i++)
            a[i] *= b[i];
        FFT(a, n, complex<double>(cos(-phi), sin(-phi)));
        BigInteger res;
        long long delta = 0;
        size_t pos = 0;
        int cntzero = 0;
        while (pos < n || delta) {
            delta += round(a[pos].real() / n);
            if (delta % MOD) {
                for (int i = 0; i < cntzero; i++)
                    res.data.push_back(0);
                cntzero = 0;
                res.data.push_back(delta % MOD);
            } else
                cntzero++;
            delta /= MOD;
            pos++;
        }
        res.sign = mulsigns(sign, second.sign);
        delete[] a;
        delete[] b;
        return res;
    }

public:
    deque<int> data;
    signs sign;

    BigInteger() : sign(signs::zero) {}

    BigInteger(long long n) {
        if (n == 0) {
            sign = signs::zero;
        } else if (n > 0) {
            sign = signs::pos;
        } else
            sign = signs::neg;
        n = std::abs(n);
        while (n > 0) {
            data.push_back(n % MOD);
            n /= MOD;
        }
    }

    BigInteger(const BigInteger &base) {
        data = base.data;
        sign = base.sign;
    }

    BigInteger &operator=(const BigInteger &second) {
        data = second.data;
        sign = second.sign;
        return *this;
    }

    BigInteger &operator++() {
        if (sign == signs::zero) {
            sign = signs::pos;
            data.push_back(1);
        } else if (sign == signs::pos) {
            add1();
        } else {
            sub1();
        }
        return *this;
    }

    BigInteger &operator--() {
        if (sign == signs::zero) {
            sign = signs::neg;
            data.push_back(1);
        } else if (sign == signs::neg) {
            add1();
        } else {
            sub1();
        }
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger res = BigInteger(*this);
        ++(*this);
        return res;
    }

    BigInteger operator--(int) {
        BigInteger res = BigInteger(*this);
        --(*this);
        return res;
    }

    BigInteger operator-() const {
        BigInteger res = BigInteger(*this);
        res.sign = mulsigns(res.sign, signs::neg);
        return res;
    }

    string toString() const {
        string s = "";
        if (sign == signs::neg)
            s += '-';
        for (int i = data.size() - 1; i >= 0; i--) {
            string s0 = "";
            int now = data[i];
            while (now) {
                s0.push_back('0' + now % 2);
                now /= 2;
            }
            std::reverse(s0.begin(), s0.end());
            if (i != static_cast<int>(data.size()) - 1) {
                while (s0.size() < DIGITS) {
                    s0 = "0" + s0;
                }
            }
            s += s0;
        }
        if (sign == signs::zero)
            s = "0";
        return s;
    }

    explicit operator bool() const {
        return sign != signs::zero;
    }

    BigInteger &operator+=(const BigInteger &second);
    BigInteger &operator-=(const BigInteger &second);
    BigInteger &operator*=(const BigInteger &second);
    BigInteger &operator/=(const BigInteger &second);
    BigInteger &operator%=(const BigInteger &second);

    friend BigInteger operator+(const BigInteger &first, const BigInteger &second);
    friend BigInteger operator-(const BigInteger &first, const BigInteger &second);
    friend std::ostream &operator<<(std::ostream &output, const BigInteger &val);
    friend std::istream &operator>>(std::istream &input, BigInteger &val);
    friend BigInteger operator*(const BigInteger &first, const BigInteger &second);
    friend BigInteger operator/(const BigInteger &first, const BigInteger &second);
    friend std::strong_ordering operator<=>(const BigInteger &first, const BigInteger &second);
    friend bool operator==(const BigInteger &first, const BigInteger &second);
    friend PoweredInteger divide(const BigInteger &first, const BigInteger &second, size_t precision);
    friend BigInteger gcd(BigInteger a, BigInteger b);
};

std::strong_ordering operator<=>(const BigInteger &first, const BigInteger &second) {
    if (first.sign == signs::pos && second.sign != signs::pos)
        return std::strong_ordering::greater;
    if (first.sign == signs::neg && second.sign != signs::neg)
        return std::strong_ordering::less;
    if (first.sign == signs::zero)
        return 0 <=> static_cast<int>(second.sign);
    const BigInteger &swapped1 = (first.sign == signs::pos) ? first : second;
    const BigInteger &swapped2 = (first.sign == signs::pos) ? second : first;
    return swapped1.cmpabs(swapped2);
}

bool operator==(const BigInteger &first, const BigInteger &second) {
    if (first.sign != second.sign)
        return false;
    if (first.data.size() != second.data.size())
        return false;
    for (int i = first.data.size() - 1; i >= 0; --i) {
        if (first.data[i] != second.data[i])
            return false;
    }
    return true;
}

bool operator!=(const BigInteger &first, const BigInteger &second) {
    return !(first == second);
}

BigInteger operator+(const BigInteger &first, const BigInteger &second) {
    return first.add(second);
}

BigInteger operator*(const BigInteger &first, const BigInteger &second) {
    return first.mul(second);
}

BigInteger operator-(const BigInteger &first, const BigInteger &second) {
    return first.add(second, signs::neg);
}
std::ostream &operator<<(std::ostream &output, const BigInteger &val) {
    output << val.toString();
    return output;
}

BigInteger p[20];

BigInteger transform(vector<char> &a) {
    int n = a.size();

    if (n == 1) {
        BigInteger res = BigInteger(a[0]);
        return res;
    }
    int logn = 0;
    while ((1 << logn) < n)
        logn++;
    vector<char> a0, a1;
    for (int i = 0; i < n / 2; i++) {
        a0.push_back(a[i]);
        a1.push_back(a[i + n / 2]);
    }
    BigInteger b0 = transform(a0);
    BigInteger b1 = transform(a1);
    return b0 + p[logn - 1] * b1;
}

int main() {
    std::ios::sync_with_stdio(0);
    cin.tie(0);
    vector<char> a;
    int base = 10;
    string s;
    cin >> s;
    for (auto z : s) {
        a.push_back(z - '0');
    }
    size_t n = 1;
    while (n < a.size())
        n *= 2;
    std::reverse(a.begin(), a.end());
    while (a.size() < n)
        a.push_back(0);
    p[0] = BigInteger(base);
    for (int i = 1; i < 20; i++) {
        p[i] = p[i - 1] * p[i - 1];
    }
    cout << transform(a);
}