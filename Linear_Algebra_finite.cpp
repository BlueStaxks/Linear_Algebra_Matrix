#include <vector>
#include <iostream>
using namespace std;
long long MOD = 10e9; //should be smaller sqrt(LLONG_MAX) for overflow

template <typename T>
inline bool operator == (vector<vector<T>> &a, vector<vector<T>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size())  return false;
    for(int i=0; i<a.size(); ++i)
        for(int j=0; j<a[0].size(); ++j)
            if(a[i][j] != b[i][j])
                return false;
    return true;
}
template <typename T>
inline bool operator != (vector<vector<T>> &a, vector<vector<T>> &b) {
    return !(a==b);
}
template <typename T>
inline bool operator == (vector<T> &a, vector<T> &b) {
    if (a.size() != b.size())  return false;
    for(int i=0; i<a.size(); ++i)
        if(a[i] != b[i])
            return false;
    return true;
}
template <typename T>
inline bool operator != (vector<T> &a, vector<T> &b) {
    return !(a==b);
}


inline vector<vector<long long>> operator * (const vector<vector<long long>> &a, const vector<vector<long long>> &b) {
    if (a.front().size() != b.size()) {
        printf("Matrix Multiplication Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(),0));
    int i, j, k;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            for (k = 0; k < b.size(); ++k) {
                R[i][j] += a[i][k] * b[k][j];
                R[i][j] %= MOD;
            }
    return R;
}
inline vector<long long> operator * (const vector<vector<long long>> &a, const vector<long long> &b) {
    if (a.front().size() != b.size()) {
        printf("Matrix Vector Multiplication Error : Matrix and Vector's size do not match\n\n");
        exit(1);
    }
    vector<long long> R(a.size());
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.size(); ++j) {
            R[i] += a[i][j] * b[j];
            R[i] %= MOD;
        }
    return R;
}
inline long long operator * (const vector<long long> &a, const vector<long long> &b) {
    if (a.size() != b.size()) {
        printf("Vector Dot Product Error : Vector size does not match\n\n");
        exit(1);
    }
    long long r=0;
    for(auto i=0; i<a.size(); ++i) {
        r += a[i]*b[i];
        r %= MOD;
    }
    return r;
}
inline vector<long long> operator * (const long long &a, const vector<long long> &b) {
    vector<long long> R(b.size());
    for(auto i=0; i<b.size(); ++i)
        R[i] = (a*b[i])%MOD;;
    return R;
}
inline vector<vector<long long>> operator + (const vector<vector<long long>> &a, const vector<vector<long long>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Addition Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(),0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = (a[i][j] + b[i][j]) % MOD;
    return R;
}
inline vector<vector<long long>> operator - (const vector<vector<long long>> &a, const vector<vector<long long>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Subtraction Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(),0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = (a[i][j] - b[i][j]) % MOD;
    return R;
}
inline vector<long long> operator + (const vector<long long> &a, const vector<long long> &b) {
    if (a.size() != b.size()) {
        printf("Vector Addition Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<long long> R(a.size(),0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = (a[i] + b[i]) % MOD;
    return R;
}
inline vector<long long> operator - (const vector<long long> &a, const vector<long long> &b) {
    if (a.size() != b.size()) {
        printf("Vector Subtraction Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<long long> R(a.size(),0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = (a[i] - b[i]) % MOD;
    return R;
}
inline long long Ratio_power(long long a, unsigned long long n) {
    long long res = 1;
    while (n) {
        if (n & 1)  res = (res * a) % MOD;
        n >>= 1;
        a = (a * a) % MOD;
    }
    return res;
}
inline long long gcd(long long a, long long b) {
    if(a<b) a^=b^=a^=b;
    long long n;
    while(b) {
        n = a % b;
        a = b;
        b = n;
    }
    return a;
}
inline vector<long long> Extended_Euclid(long long a, long long b) {
    long long q,r1=a,r2=b,r=1,s1=1,s2=0,s,t1=0,t2=1,t;
    while(r) {
        q=r1/r2;    r=r1%r2;
        s=(s1-q*s2);
        t=(t1-q*t2);
        r1=r2;  r2=r;   s1=s2;  s2=s;   t1=t2;  t2=t;
    }
    return {r1,s1,t1};
}
inline long long int_inverse(long long a) {
    if(!a) {
        printf("Integer Inverse Error : 0 has no inverse.\n\n");
        exit(1);
    }
    long long q,r1=MOD,r2=a,r=1,t1=0,t2=1,t;
    while(r) {
        q=r1/r2;    r=r1%r2;
        t=(t1-q*t2) % MOD;
        if(t<0) t+=MOD;
        r1=r2;  r2=r;   t1=t2;  t2=t;
    }
    if(r1!=1) {
        printf("Integer Inverse Error : %lld does not have inverse in Z(%lld)\n\n",a,MOD);
        exit(1);
    }
    return t1;
}
int main()
{
    MOD = 100;
    long long k = int_inverse(23);
    printf("%lld\n",k);
    vector<long long> v = Extended_Euclid(161, 28);
    for(auto i: v)
        printf("%lld ",i);
    return 0;
}
