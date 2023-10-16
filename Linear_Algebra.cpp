#include <vector>
#include <math.h>
#include <iostream>
using namespace std;
inline unsigned long long gcd(unsigned long long a, unsigned long long b) {
    if(a<b) a^=b^=a^=b;
    unsigned long long n;
    while(b) {
        n = a % b;
        a = b;
        b = n;
    }
    return a;
}
typedef struct Ratio {
    bool sign; //true N   false P
    unsigned long long a,b;
    
    Ratio() {a=0;b=1;sign=false;}
    Ratio(long long k) {
        if(k<0) {
            sign=true;
            a=-k;
        }
        else {
            sign=false;
            a=k;
        }
        b=1;
    }
    Ratio(long long k1, long long k2) {
        if(!k2) {
            printf("Ratio constructor Error : denominator is zero.\n\n");
            exit(1);
        }
        if(k1<0 && k2<0) {
            sign=false;
            a=-k1;  b=-k2;
        }
        else if(k1>=0 && k2>=0) {
            sign=false;
            a=k1;   b=k2;
        }
        else if(k1<0 && k2>=0) {
            sign=true;
            a=-k1;  b=k2;
        }
        else {
            sign=true;
            a=k1;   b=-k2;
        }
    }
    Ratio(bool S, unsigned long long k1, unsigned long long k2) {
        if(!k2) {
            printf("Ratio constructor Error : denominator is zero.\n\n");
            exit(1);
        }
        sign=S; a=k1;   b=k2;
    }
    inline void normalize()
    {
        if(!a) {
            sign = false;
            b=1;
            return;
        }
        unsigned long long d = gcd(a,b);
        a/=d;   b/=d;
    }
    inline bool operator < (const Ratio &K) {
        if(sign && !K.sign) return true;
        else if(!sign && K.sign)    return false;
        unsigned long long d=b*K.b/gcd(b,K.b);
        bool r = a*d/b < K.a*d/K.b;
        if(sign)    return !r;
        else    return r;
    }
    inline bool operator <= (const Ratio &K) {
        if(*this==K)    return true;
        if(*this < K)   return true;
        return false;
    }
    inline bool operator > (const Ratio &K) {
        if(sign && !K.sign) return false;
        else if(!sign && K.sign)    return true;
        unsigned long long d=b*K.b/gcd(b,K.b);
        bool r = a*d/b > K.a*d/K.b;
        if(sign)    return !r;
        else    return r;
    }
    inline bool operator >= (const Ratio &K) {
        if(*this==K)    return true;
        if(*this > K)   return true;
        return false;
    }
    inline bool operator == (const Ratio &K) {
        return a==K.a && b==K.b && sign==K.sign;
    }
    inline bool operator != (const Ratio &K) {
        return a!=K.a || b!=K.b || sign!=K.sign;
    }
    inline void Ratio_print(bool F) {
        if(sign)    printf("-");
        if(F)   printf("%Lf\t",(long double)a/(long double)b);
        else    printf("%llu / %llu\t\t",a,b);
    }
} Ratio;
inline Ratio operator * (Ratio A, Ratio B) {
    if(!A.b || !B.b) {
        printf("Ratio Error : divide by 0\n\n");
        exit(1);
    }
    if(!A.a || !B.a)    return {false, 0,1};
    unsigned long long d1 = gcd(A.a,B.b);
    A.a/=d1;  B.b/=d1;
    unsigned long long d2 = gcd(B.a,A.b);
    B.a/=d2;  A.b/=d2;
    return {static_cast<bool>(A.sign^B.sign), A.a*B.a, A.b*B.b};
}
inline Ratio operator * (Ratio A, long long B) {
    if(!A.b) {
        printf("Ratio Error : divide by 0\n\n");
        exit(1);
    }
    if(!A.a)    return {false, 0,1};
    unsigned long long d2 = gcd(abs(B),A.b);
    B/=d2;  A.b/=d2;
    unsigned long long t = abs(B);
    return {static_cast<bool>(A.sign^(B<0)), static_cast<long long>(A.a*t)};
}
inline Ratio operator / (Ratio A, Ratio B) {
    B.a^=B.b^=B.a^=B.b;
    return A*B;
}
inline Ratio operator + (Ratio A, Ratio B) {
    if(!A.a)    return B;
    if(!B.a)    return A;
    unsigned long long d = gcd(A.b,B.b);
    unsigned long long t = A.b*B.b/d;
    unsigned long long c1 = B.b/d, c2 = A.b/d;
    A.a*=c1;    B.a*=c2;    A.b=B.b=t;
    if(!A.sign && !B.sign) {
        Ratio r = {false, A.a+B.a, t};
        r.normalize();
        return r;
    }
    if(!A.sign && B.sign) {
        if(A.a >= B.a) {
            Ratio r = {false, A.a-B.a, t};
            r.normalize();
            return r;
        }
        else {
            Ratio r = {true, B.a-A.a, t};
            r.normalize();
            return r;
        }
    }
    if(A.sign && !B.sign) {
        if(A.a <= B.a) {
            Ratio r = {false, B.a-A.a, t};
            r.normalize();
            return r;
        }
        else {
            Ratio r = {true, A.a-B.a, t};
            r.normalize();
            return r;
        }
    }
    else {
        Ratio r = {true, A.a+B.a, t};
        r.normalize();
        return r;
    }
}
inline Ratio operator - (Ratio A, Ratio B) {
    B.sign = !B.sign;
    return A+B;
} // -------------------------------------------------------------------------- Ratio define
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
template <>
inline bool operator == (vector<vector<long double>> &a, vector<vector<long double>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size())  return false;
    for(int i=0; i<a.size(); ++i)
        for(int j=0; j<a[0].size(); ++j)
            if(abs(a[i][j] - b[i][j]) > 0.00000001)
                return false;
    return true;
}
inline bool operator == (vector<vector<Ratio>> &a, vector<vector<long double>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size())  return false;
    for(int i=0; i<a.size(); ++i)
        for(int j=0; j<a[0].size(); ++j) {
            long double ld = a[i][j].a / (long double)a[i][j].b;
            if(a[i][j].sign)    ld=-ld;
            if(abs(ld - b[i][j]) > 0.000001)
                return false;
        }
    return true;
}
inline bool operator != (vector<vector<Ratio>> &a, vector<vector<long double>> &b) {
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
template <>
inline bool operator == (vector<long double> &a, vector<long double> &b) {
    if (a.size() != b.size())  return false;
    for(int i=0; i<a.size(); ++i)
        if(abs(a[i] - b[i]) > 0.000001)
            return false;
    return true;
}
template <typename T>
inline vector<vector<T>> operator * (const vector<vector<T>> &a, const vector<vector<T>> &b) {
    if (a.front().size() != b.size()) {
        printf("Matrix Multiplication Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<T>> R(a.size(), vector<T>(b.front().size(),0));
    int i, j, k;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            for (k = 0; k < b.size(); ++k)
                R[i][j] = R[i][j] + (a[i][k] * b[k][j]);
    return R;
}
template <typename T>
inline vector<T> operator * (const vector<vector<T>> &a, const vector<T> &b) {
    if (a.front().size() != b.size()) {
        printf("Matrix Vector Multiplication Error : Matrix and Vector's size do not match\n\n");
        exit(1);
    }
    vector<T> R(a.size());
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.size(); ++j)
            R[i] = R[i] + (a[i][j] * b[j]);
    return R;
}
template <typename T>
inline T operator * (const vector<T> &a, const vector<T> &b) {
    if (a.size() != b.size()) {
        printf("Vector Dot Product Error : Vector size does not match\n\n");
        exit(1);
    }
    T r=0;
    for(int i=0; i<a.size(); ++i)
        r = r + a[i]*b[i];
    return r;
}
template <typename T>
inline vector<T> operator * (const T &a, const vector<T> &b) {
    vector<T> R(b.size());
    for(int i=0; i<b.size(); ++i)
        R[i] = a*b[i];
    return R;
}
template <typename T>
inline vector<vector<T>> operator + (const vector<vector<T>> &a, const vector<vector<T>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Addition Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<T>> R(a.size(), vector<T>(b.front().size(),0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = a[i][j] + b[i][j];
    return R;
}
template <typename T>
inline vector<vector<T>> operator - (const vector<vector<T>> &a, const vector<vector<T>> &b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Subtraction Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<T>> R(a.size(), vector<T>(b.front().size(),0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = a[i][j] - b[i][j];
    return R;
}
template <typename T>
inline vector<T> operator + (const vector<T> &a, const vector<T> &b) {
    if (a.size() != b.size()) {
        printf("Vector Addition Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<T> R(a.size(),0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = a[i] + b[i];
    return R;
}
template <typename T>
inline vector<T> operator - (const vector<T> &a, const vector<T> &b) {
    if (a.size() != b.size()) {
        printf("Vector Subtraction Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<T> R(a.size(),0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = a[i] - b[i];
    return R;
}
template <typename T>
inline vector<T> Vector_Normalize(vector<T> v) {
    int i;
    T s=0;
    for(i=0; i<v.size(); ++i)   s+=v[i]*v[i];
    s = sqrt(s);
    for(i=0; i<v.size(); ++i)   v[i]/=s;
    return v;
}
template <typename T>
inline vector<T> Vector_Denormalize(vector<T> v) {
    int i;
    T min = __LDBL_MAX__;
    vector<T> r(v.size(),0);
    for(i=0; i<v.size(); ++i) {
        if(v[i]==0) continue;
        r[i]=round(v[i]*134217728);
        if(abs(r[i]) < min) min = abs(r[i]);
    }
    for(i=0; i<r.size(); ++i)   r[i]/=min;
    return r;
}
template <typename T>
inline vector<vector<T>> matrix_transpose(const vector<vector<T>> &a) {
    vector<vector<T>> R(a.front().size(), vector<T>(a.size()));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < a.front().size(); ++j)
            R[j][i] = a[i][j];
    return R;
}
inline vector<vector<long double>> RatioMat_to_LDMat(const vector<vector<Ratio>>& A) {
    vector<vector<long double>> R(A.size(), vector<long double>(A[0].size()));
    for(int i=0; i<A.size(); ++i)
        for(int j=0; j<A[i].size(); ++j) {
            R[i][j] = A[i][j].a / (long double)A[i][j].b;
            if(A[i][j].sign)    R[i][j]=-R[i][j];
        }
    return R;
}
template <typename T>
inline vector<vector<T>> matrix_power(vector<vector<T>> a, unsigned long long n) {
    auto size = a.size();
    vector<vector<T>> res(size, vector<T>(size, 0));
    for (int i = 0; i < size; i++)  res[i][i] = 1;
    while (n) {
        if (n & 1)  res = res * a;
        n >>= 1;
        a = a * a;
    }
    return res;
}
template <typename T>
inline vector<vector<T>> matrix_row_multiply(vector<vector<T>> A, const vector<T>& v) {
    if(A.size() != v.size())    exit(1);
    for(int i=0; i<A.size(); ++i)
        for(int j=0; j<A[i].size(); ++j)
            A[i][j] = A[i][j] * v[i];
    return A;
}
template <typename T>
inline T Ratio_power(T a, unsigned long long n) {
    T res = 1;
    while (n) {
        if (n & 1)  res = res * a;
        n >>= 1;
        a = a * a;
    }
    return res;
}
inline void matrix_print(vector<vector<Ratio>> a, bool F) {
    if(a.empty())   return;   if(a[0].empty())    return;
    for(int i=0; i<a.size(); ++i)
    {
        for(int j=0; j<a[0].size(); ++j)
            a[i][j].Ratio_print(F);
        printf("\n");
    }
    printf("\n\n");
}
inline void matrix_print(vector<vector<long double>> a) {
    if(a.empty())   return;   if(a[0].empty())    return;
    for(int i=0; i<a.size(); ++i)
    {
        for(int j=0; j<a[0].size(); ++j)
            printf("%Lf\t",a[i][j]);
        printf("\n");
    }
    printf("\n\n");
}
inline void vector_print(vector<Ratio>& a, bool F) {
    for(int i=0; i<a.size(); ++i)
        a[i].Ratio_print(F);
    printf("\n\n\n");
}
inline void vector_print(vector<long double>& a) {
    for(int i=0; i<a.size(); ++i)
        printf("%Lf\t",a[i]);
    printf("\n\n\n");
}
template <typename T>
inline vector<vector<T>> matrix_inverse(vector<vector<T>> A) {
    if (A.size() != A.front().size()) {
        printf("Matrix Inversion Error : Matrix is not square\n\n");
        exit(1);
    }
    auto n = A.size();
    long long i, j, k;
    vector<vector<T>> I(n, vector<T>(n, 0));
    for (i = 0; i < n; ++i)  I[i][i] = 1;
    for (i = 1; i < n; ++i) {
        if (A[i - 1][i - 1]==0)
        {
            bool P=true;
            for(j=i; j<n; ++j)
                if(A[j][i-1]!=0)
                {
                    vector<T> t=A[j];   A[j]=A[i-1];    A[i-1]=t;
                    t=I[j]; I[j]=I[i-1];    I[i-1]=t;
                    P=false;
                    break;
                }
            if(P) {
                printf("Matrix Inversion Error : Matrix is singular\n\n");
                exit(1);
            }
        }
        for (j = i; j < n; ++j) {
            T mul = A[j][i - 1] / A[i - 1][i - 1];
            for (k = 0; k < n; ++k)
            {
                A[j][k] = A[j][k] - (A[i - 1][k] * mul);
                I[j][k] = I[j][k] - (I[i - 1][k] * mul);
            }
        }
    }
    for (i = n - 2; i >= 0; --i) {
        if (A[i + 1][i + 1]==0) {
            printf("Matrix Inversion Error : Matrix is singlular\n\n");
            exit(1);
        }
        for (j = i; j >= 0; --j) {
            T mul = A[j][i + 1] / A[i + 1][i + 1];
            for (k = 0; k < n; ++k)
            {
                A[j][k] = A[j][k] - (A[i + 1][k] * mul);
                I[j][k] = I[j][k] - (I[i + 1][k] * mul);
            }
        }
    }
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            I[i][j] = I[i][j] / A[i][i];
    return I;
}
template <typename T>
inline T matrix_determinant(vector<vector<T>> a) {
    if (a.size() != a.front().size()) {
        printf("Matrix determinant Error : Matrix is not square\n\n");
        exit(1);
    }
    int tr = 1;
    auto n = a.size();
    long long i, j, k;
    for (i = 1; i < n; ++i) {
        if (a[i - 1][i - 1]==0)
        {
            bool P=true;
            for(j=i; j<n; ++j)
                if(a[j][i-1]!=0)
                {
                    vector<T> t=a[j];   a[j]=a[i-1];    a[i-1]=t;
                    tr*=-1;
                    P=false;
                    break;
                }
            if(P)   return 0;
        }
        for (j = i; j < n; ++j) {
            T mul = a[j][i - 1] / a[i - 1][i - 1];
            if(!mul.a)  continue;
            for (k = 0; k < n; ++k)
                a[j][k] = a[j][k] - (a[i - 1][k] * mul);
        }
    }
    T r = tr;
    for(i=0; i<n; ++i)  r = r * a[i][i];
    return r;
}
template <typename T>
inline void QR_decomposition(vector<vector<T>>& A, vector<vector<T>>& Q, vector<vector<T>>& R) {
    R.resize(A[0].size(), vector<T>(A[0].size(),0));    R[0][0]=1;
    vector<vector<T>> X = matrix_transpose(A);
    vector<vector<T>> V = X;
    vector<T> DP(2,0);
    vector<T> N(A[0].size());
    int i,j;
    DP[1] = V[0]*V[0];
    for(i=1; i<A[0].size(); ++i)
    {
        for(j=1; j<=i; ++j)
        {
            T c = (V[j-1]*X[i])/DP[j];
            R[j-1][i] = c;
            V[i] = V[i] - (c*V[j-1]);
        }
        DP.push_back(V[i]*V[i]);
        R[j-1][i] = (X[i]*V[i])/DP.back();
        if(DP.back()==0) {
            printf("QR decomposition Error : Matrix is singular\n\n");
            exit(1);
        }
    }
    for(i=0; i<V.size(); ++i)
    {
        T r = 0;
        for(j=0; j<V[0].size(); ++j)    r += V[i][j] * V[i][j];
        r = sqrt(r);
        for(j=0; j<V[0].size(); ++j)    V[i][j] /= r;
        N[i]=r;
    }
    Q = matrix_transpose(V);
    R = matrix_row_multiply(R, N);
}
template <>
inline void QR_decomposition(vector<vector<Ratio>>& A, vector<vector<Ratio>>& Q, vector<vector<Ratio>>& R) {
    R.resize(A[0].size(), vector<Ratio>(A[0].size(),0));    R[0][0]=1;
    vector<vector<Ratio>> X = matrix_transpose(A);
    vector<vector<Ratio>> V = X;
    vector<Ratio> DP(2,0);
    int i,j;
    DP[1] = V[0]*V[0];
    for(i=1; i<A[0].size(); ++i)
    {
        for(j=1; j<=i; ++j)
        {
            Ratio c = (V[j-1]*X[i])/DP[j];
            R[j-1][i] = c;
            V[i] = V[i] - (c*V[j-1]);
        }
        DP.push_back(V[i]*V[i]);
        R[j-1][i] = (X[i]*V[i])/DP.back();
        if(DP.back()==0) {
            printf("QR decomposition Error : Matrix is singular\n\n");
            exit(1);
        }
    }
    Q = matrix_transpose(V);
}
template <typename T>
inline void LU_decomposition(vector<vector<T>> A, vector<vector<T>>& L, vector<vector<T>>& U) {
    auto m = A.size(), n = A[0].size();
    L.resize(m,vector<T>(m,0));
    auto el = n>m?m:n;
    int i,j,k,p=0;
    for(i=0; i<m; ++i)  L[i][i]=1;
    for(i=1; i<=el; ++i)
    {
        //matrix_print(A, 0); printf("\n\n");
        if(A[i-1-p][i-1]==0)
        {
            for(j=i-p; j<m; ++j)
                if(A[j][i-1]!=0) {
                   printf("LU decompotision Error : Zero pivot\n\n");
                   exit(1);
                }
            p++;
            continue;
        }
        for (j = i - p; j < m; ++j) {
            T mul = A[j][i - 1] / A[i - 1 - p][i - 1];
            L[j][i-1-p] = mul;
            if(!mul.a)  continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = A[j][k] - (A[i - 1 - p][k] * mul);
        }
    }
    for(i=(int)n-p; i<m; ++i)  L[i][i]=1;
    U.resize(m,vector<T>(n));
    for(i=0; i<m; ++i)
        for(j=0; j<n; ++j)
            U[i][j]=A[i][j];
}
template <typename T>
inline vector<vector<T>> clean_eigenvector(vector<vector<T>>& S) {
    vector<vector<T>> TR = matrix_transpose(S), R(S.size(), vector<T>(S[0].size()));
    for(int i=0; i<TR.size(); ++i)  R[i]=Vector_Denormalize(TR[i]);
    return matrix_transpose(R);
}
template <typename T>
inline void Eigen_Approx(vector<vector<T>>& A, vector<vector<T>>& e_val, vector<vector<T>>& e_vec, int n) {
    if(A.size() != A[0].size()) {
        printf("Eigen Approximation Error : Matrix is not square\n\n");
        exit(1);
    }
    for(int i=0; i<A.size(); ++i)
        for(int j=0; j<i; ++j)
            if(A[i][j] != A[j][i]) {
                printf("Eigen Approximation Error : Matrix is not Symmetric\n\n");
                exit(1);
            }
    vector<vector<T>> Q, R;
    e_val.resize(A.size(), vector<T>(A[0].size()));
    e_vec.resize(A.size(), vector<T>(A[0].size(),0));
    for(int i=0; i<A.size(); ++i)
    {
        for(int j=0; j<A[0].size(); ++j)
            e_val[i][j]=A[i][j];
        e_vec[i][i]=1;
    }
    for(int i=0; i<n; ++i) {
        QR_decomposition(e_val, Q, R);
        e_vec = e_vec * Q;
        e_val = R * Q;
    }
    QR_decomposition(e_val, Q, R);
    e_vec = e_vec * Q;
    vector<vector<T>> K = e_vec * e_val * matrix_transpose(e_vec);
    if(A!=K) {
        printf("Eigen Approximation Failed : Maybe not enough n\n\n");
        exit(1);
    }
    return;
}
template <>
inline void Eigen_Approx(vector<vector<Ratio>>& A, vector<vector<Ratio>>& e_val, vector<vector<Ratio>>& e_vec, int n) {
    printf("Eigen_Approx is not possible with Ratio type Matrix\n\n");
    exit(1);
}
template <typename T>
inline vector<vector<T>> Null_Space(vector<vector<T>> A) {
    int m = (int)A.size(), n = (int)A[0].size();
    vector<int> piv;
    int i,j,k,p=0,rank=0;
    for(i=1; i-1<n && i-1-p<m; ++i)
    {
        if(A[i-1-p][i-1]==0)
        {
            bool P=true;
            for(j=i-p; j<m; ++j)
                if(A[j][i-1]!=0) {
                    vector<T> temp = A[i-1-p];
                    A[i-1-p] = A[j];
                    A[j] = temp;
                    i--;
                    P=false;
                    break;
                }
            if(!P)  continue;
            p++;
            continue;
        }
        piv.push_back(i-1);
        rank++;
        T temp = A[i-1-p][i-1];
        A[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)    A[i-1-p][j] = A[i-1-p][j] / temp;
        for (j = i - p; j < m; ++j) {
            T mul = A[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = A[j][k] - (A[i - 1 - p][k] * mul);
        }
    }
    for(i=(int)piv.size()-1; i>0; --i) //upper elimination
        for(j=i-1; j>=0; --j)
        {
            T mul = A[j][piv[i]];
            for(k=piv[i]; k<n; ++k)
                A[j][k] = A[j][k] - (A[i][k] * mul);
        }
    for(i=m-1; i>=0; --i) //zero row pop
    {
        bool P=false;
        for(j=0; j<n; ++j)
            if(A[i][j]!=0) {
                P=true;
                break;
            }
        if(P)   break;
        A.pop_back();
    }
    if(A.size() == A[0].size()) {
        printf("Null_Space calculation Alert : There is no Special Solutions\n\n");
        exit(1);
    }
    vector<vector<T>> TR = matrix_transpose(A), F(n-rank);
    vector<pair<int,int>> exc;
    T mo = -1;
    for(i=0; i<rank; ++i) {
        if(piv[i]!=i) {
            vector<T> te = TR[i];
            TR[i] = TR[piv[i]];
            TR[piv[i]] = te;
            exc.push_back({i,piv[i]});
            piv[i]=i; //data loss
        }
    }
    for(p=i; i<n; ++i) F[i-p] = mo * TR[i];
    vector<vector<T>> N = matrix_transpose(F);
    for(i=0; i<F.size(); ++i) {
        vector<T> te(F.size(),0);
        te[i]=1;
        N.push_back(te);
    }
    for(i=(int)exc.size()-1; i>=0; --i) {
        vector<T> te = N[exc[i].first];
        N[exc[i].first] = N[exc[i].second];
        N[exc[i].second] = te;
    }
    return N;
}
template <typename T>
inline vector<T> Ax_b(vector<vector<T>>& A, vector<T> b) {
    if(A.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }
    int m = (int)A.size(), n = (int)A[0].size();
    n++;
    vector<vector<T>> R(m,vector<T>(n,0));
    vector<int> piv;
    int i,j,k,p=0,elimi=0;
    for(i=0; i<m; ++i)
    {
        for(j=0; j<n-1; ++j)
            R[i][j]=A[i][j];
        R[i][j]=b[i];
    }
    for(i=1; i<n && i-1-p<m; ++i)
    {
        if(R[i-1-p][i-1]==0)
        {
            bool P=true;
            for(j=i-p; j<m; ++j)
                if(R[j][i-1]!=0) {
                    vector<T> temp = R[i-1-p];
                    R[i-1-p] = R[j];
                    R[j] = temp;
                    i--;
                    P=false;
                    break;
                }
            if(!P)  continue;
            p++;
            continue;
        }
        piv.push_back(i-1);
        elimi++;
        T temp = R[i-1-p][i-1];
        R[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)    R[i-1-p][j] = R[i-1-p][j] / temp;
        for (j = i - p; j < m; ++j) {
            T mul = R[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                R[j][k] = R[j][k] - (R[i - 1 - p][k] * mul);
        }
    }
    for(i=(int)piv.size()-1; i>0; --i) //upper elimination
        for(j=i-1; j>=0; --j)
        {
            T mul = R[j][piv[i]];
            for(k=piv[i]; k<n; ++k)
                R[j][k] = R[j][k] - (R[i][k] * mul);
        }
    for(i=m-1; i>=0; --i) //zero row solvablity
    {
        bool P=false;
        for(j=0; j<n-1; ++j)
            if(R[i][j]!=0) {
                P=true;
                break;
            }
        if(P)   break;
        if(R[i][n-1]!=0) {
            printf("Ax=b calculation Error : This System is Not Solvable\n\n");
            exit(1);
        }
        R.pop_back();
        b.pop_back();
    }
    if(R.size() == R[0].size()-1) {
        vector<T> r(R.size());
        for(i=0; i<r.size(); ++i)   r[i]=R[i][n-1];
        return r;
    }
    vector<T> r(n-1,0);
    for(i=0, p=0; i<n-1; ++i) {
        if(piv[p]==i) {
            r[i]=R[p][n-1];
            p++;
        }
    }
    return r;
}
template <typename T>
inline vector<vector<T>> matrix_full_row_rank(vector<vector<T>> A) {
    auto m = A.size(), n = A[0].size();
    vector<vector<T>> A2 = A;
    vector<int> exc(m);
    vector<bool> Ap(m,false);
    int i,j,k,p=0;
    for(i=0; i<m; ++i)  exc[i]=i;
    for(i=1; i-1<n && i-1-p<m; ++i)
    {
        if(A[i-1-p][i-1]==0)
        {
            bool P=true;
            for(j=i-p; j<m; ++j)
                if(A[j][i-1]!=0) {
                    exc[i-1-p] ^= exc[j] ^= exc[i-1-p] ^= exc[j];
                    vector<T> temp = A[i-1-p];
                    A[i-1-p] = A[j];
                    A[j] = temp;
                    i--;
                    P=false;
                    break;
                }
            if(!P)  continue;
            p++;
            continue;
        }
        T temp = A[i-1-p][i-1];
        A[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)    A[i-1-p][j] = A[i-1-p][j] / temp;
        for (j = i - p; j < m; ++j) {
            T mul = A[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = A[j][k] - (A[i - 1 - p][k] * mul);
        }
    }
    for(i=m-1; i>=0; --i) //zero row pop
    {
        bool P=false;
        for(j=0; j<n; ++j)
            if(A[i][j]!=0) {
                P=true;
                break;
            }
        if(P)   break;
        //A.pop_back();
        Ap[exc[i]]=true;
    }
    vector<vector<T>> FRA;
    for(i=0; i<m; ++i)
        if(!Ap[i])
            FRA.push_back(A2[i]);
    return FRA;
}
template <typename T>
inline vector<vector<T>> change_of_basis_P(vector<vector<T>> B, vector<vector<T>> C) { //B to C
    if(matrix_determinant(B)==0 || matrix_determinant(C)==0) {
        printf("Change of Basis P Error : Vectors are dependent.\n\n\n");
        exit(1);
    }
    auto n = B.size();
    vector<vector<T>> R, INV = matrix_inverse(C), BT = matrix_transpose(B);
    for(int i=0; i<n; ++i)
        R.push_back(INV * BT[i]);
    return matrix_transpose(R);
}
int main()
{
    vector<vector<Ratio>> A = {
//        {1,3,0,2,-1},
//        {0,0,1,4,-3},
//        {1,3,1,6,-4} 
        
//        {1,2,2,2},
//        {2,4,6,8},
//        {3,6,8,10}
        
//        {1,2,2,2,2,8},
//        {2,4,6,82,2,4},
//        {330,6,8,30,9991,9}
        
        {1,2,3},
        {2,5,9},
        {3,9,88}
        
        
//        {1,1},
//        {3,3}
    }, L,U,Q,R;
//    vector<Ratio> b = {1,2,30};
//    QR_decomposition(A, Q, R);
//    matrix_print(Q,0);
//    matrix_print(matrix_transpose(Q) * Q,0);
    
//    vector<vector<long double>> eval, evec, LD_A;
//    LD_A = RatioMat_to_LDMat(A);
//    Eigen_Approx(LD_A, eval, evec, 10000);
//    matrix_print(eval);
//    
//    vector<vector<long double>> vec2 = matrix_transpose(evec);
////    for(int i=0; i<vec2.size(); ++i) {
////        vector<long double> vt = Vector_Denormalize(vec2[i]);
////        vector_print(vt);
////    }
//    matrix_print(evec);
//
//    matrix_print(evec * eval * matrix_transpose(evec));
//
//    vector<vector<long double>> T = matrix_transpose(evec);
//    for(int i=0; i<T.size() - 1; ++i)
//        for(int j=i+1; j<T.size(); ++j)
//            printf("%Lf\n",T[i]*T[j]);
    
//    vector<Ratio> r = Ax_b(A, b);
//    vector_print(r, 0);
//    vector<Ratio> b2 = A*r;
//    vector_print(b2, 0);
//    
//    vector<vector<Ratio>> NS = Null_Space(A);
//    matrix_print(NS, 0);
//    matrix_print(A * NS, 0);
    
//    vector<vector<Ratio>> A2 = matrix_full_row_rank(A);
//    vector<vector<Ratio>> A3 = matrix_full_row_rank(matrix_transpose(A2));
//    matrix_print(A3, 0);
    
    vector<vector<Ratio>> B = {
        {1,1,-1},
        {0,1,0},
        {0,0,1}
    };
    vector<vector<Ratio>> C = {
        {2,-1,1},
        {-1,3,0},
        {1,1,0}
    };
    vector<vector<Ratio>> R1 = change_of_basis_P(B, C);
    vector<vector<Ratio>> R2 = change_of_basis_P(C, B);
    vector<vector<Ratio>> R3 = R1 * R2;
    matrix_print(R1, 0);
    matrix_print(R2, 0);
    matrix_print(R3, 0);

    return 0;
}
