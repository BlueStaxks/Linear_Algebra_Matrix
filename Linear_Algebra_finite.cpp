#include <iostream>
#include <fstream>
#include <time.h>

#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;
long long MOD = 10e9; //should be inside of INT32 range. It's GF so MOD must be a prime.
long long primitive;
vector<int> int_inverse;
vector<int> seeds; // primitive ^ seeds[i] = i
vector<long long> MOD_decompose; 
vector<long long> MOD_divisors;
vector<vector<int>> ones_roots; // 1^(1/i) = ones_roots.front() ~ back()
inline long long inverse(long long a) {
    if(!a) {
        printf("Integer Inverse Error : 0 has no inverse.\n\n");
        exit(1);
    }
    if(int_inverse[a])  return int_inverse[a];
    long long q,r1=MOD,r2=a,r=1,t1=0,t2=1,t;
    while(r) {
        q=r1/r2;    r=r1%r2;
        t=(t1-q*t2) % MOD;
        if(t<0) t+=MOD;
        r1=r2;  r2=r;   t1=t2;  t2=t;
    }
//    if(r1!=1) {
//        printf("Integer Inverse Error : %lld does not have inverse in Z(%lld)\n\n",a,MOD);
//        exit(1);
//    }
    int_inverse[a]=(int)t1;
    return t1;
}
inline long long power(long long a, long long n) {
    long long y = 1;
    while (n > 0) {
        if (n & 1)
            y = (a * y) % MOD;
        a = (a * a) % MOD;
        n>>= 1;
    }
    return y;
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
inline vector<long long> decompose(long long a) {
    vector<long long> r;
    while(!(a&1)) {
        r.push_back(2);
        a>>=1;
    }
    for(long long i=3; i<=sqrt(a); i+=2)
        if(a%i==0) {
            r.push_back(i);
            a/=i;
        }
    r.push_back(a);
    return r;
}
inline vector<long long> divisor(long long a) {
    vector<long long> r;
    long long sq = sqrt(a);
    for(long long i = 1; i <= sq; ++i)
        if (a % i == 0) {
            r.push_back(i);
            if (i*i!=a) r.push_back(a / i);
        }
    sort(r.begin(), r.end());
    return r;
}
void Initiation() {
    int_inverse.resize(MOD,0);
    MOD_decompose = decompose(MOD-1);
    MOD_divisors = divisor(MOD-1);
    primitive=1;
    for(int i=2; i<MOD; ++i) {
        bool P=true;
        for(int j=1; j<MOD_divisors.size()-1; ++j)
            if(power(i,MOD_divisors[j])==1) {
                P=false;
                break;
            }
        if(P) {
            primitive=i; //find smallest primitive root
            break;
        }
    }
    seeds.resize(MOD);
    ones_roots.resize(MOD);
    for(long long i=1,t=primitive; i<MOD; ++i, t=t*primitive%MOD) {
        seeds[t]=(int)i;  //primitive ^ seeds[i] = i
        for(int j=0; j<MOD_divisors.size()-1; ++j) // 1^(1/i) = ones_roots.front() ~ back()
            if(power(i,MOD_divisors[j])==1)
                ones_roots[MOD_divisors[j]].push_back((int)i); //calculate order of all numbers
    }
}

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
inline void matrix_print(vector<vector<long long>> a) {
    if(a.empty())   return;   if(a[0].empty())    return;
    for(int i=0; i<a.size(); ++i)
    {
        for(int j=0; j<a[0].size(); ++j)
            printf("%lld\t",a[i][j]);
        printf("\n");
    }
    printf("\n\n");
}
inline void vector_print(vector<long long> a) {
    for(int i=0; i<a.size(); ++i)
        printf("%lld\t",a[i]);
    printf("\n\n\n");
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
            for (k = 0; k < b.size(); ++k)
                R[i][j] = (R[i][j] + a[i][k] * b[k][j]) % MOD;
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
    for(auto i=0; i<a.size(); ++i)
        r = (r + a[i]*b[i]) % MOD;
    return r;
}
inline vector<long long> operator * (const long long &a, const vector<long long> &b) {
    vector<long long> R(b.size());
    for(auto i=0; i<b.size(); ++i)
        R[i] = (a*b[i])%MOD;;
    return R;
}
inline vector<vector<long long>> operator * (const long long &a, const vector<vector<long long>> &b) {
    vector<vector<long long>> R(b.size(), vector<long long>(b.front().size()));
    for(auto i=0; i<b.size(); ++i)
        for(auto j=0; j<b.front().size(); ++j)
            R[i][j] = (a*b[i][j])%MOD;;
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
            R[i][j] = (a[i][j] - b[i][j] + MOD) % MOD;
    return R;
}
inline vector<vector<long long>> operator | (const vector<vector<long long>> &a, const vector<vector<long long>> &b) { //diagonal expansion
    if(a.empty())   return b;
    if(b.empty())   return a;
    vector<vector<long long>> R(a.size()+b.size(), vector<long long>(a.front().size()+b.front().size(),0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < a.front().size(); ++j)
            R[i][j] = a[i][j];
    for(i=(int)a.size(); i<R.size(); ++i)
        for(j=(int)a.front().size(); j<R.front().size(); ++j)
            R[i][j] = b[i-a.size()][j-a.front().size()];
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
inline vector<vector<long long>> I_n(int n) {
    vector<vector<long long>> I(n, vector<long long>(n,0));
    for (int i = 0; i < n; ++i)  I[i][i] = 1;
    return I;
}
inline vector<vector<long long>> I_n(int n, long long a) {
    vector<vector<long long>> I(n, vector<long long>(n,0));
    for (int i = 0; i < n; ++i)  I[i][i] = a;
    return I;
}
inline vector<vector<long long>> matrix_transpose(const vector<vector<long long>> &a) {
    vector<vector<long long>> R(a.front().size(), vector<long long>(a.size()));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < a.front().size(); ++j)
            R[j][i] = a[i][j];
    return R;
}
inline vector<vector<long long>> matrix_power(vector<vector<long long>> a, unsigned long long n) {
    auto size = a.size();
    vector<vector<long long>> res(size, vector<long long>(size, 0));
    for (int i = 0; i < size; i++)  res[i][i] = 1;
    while (n) {
        if (n & 1)  res = res * a;
        n >>= 1;
        a = a * a;
    }
    return res;
}
inline void matrix_chop(vector<vector<vector<long long>>>& M, vector<vector<long long>>& F, vector<int> list) {
    vector<vector<long long>> N;
    for(int i=0,p=0; i<list.size(); ++i) {
        M.push_back(N);
        M.back().resize(list[i], vector<long long>(list[i]));
        for(int j=0; j<list[i]; ++j)
            for(int k=0; k<list[i]; ++k)
                M.back()[j][k]=F[j+p][k+p];
        p+=list[i];
    }
}
inline long long matrix_rank(vector<vector<long long>> A) {
    int m = (int)A.size(), n = (int)A[0].size();
    int i,j,k,l,p=0;
    long long rank=0;
    for(i=1; i-1<n && i-1-p<m; ++i)
    {
        if(A[i-1-p][i-1]==0) //pivot is zero
        {
            bool P=true;
            for(j=i-p; j<m; ++j) //row exchange is allowed
                if(A[j][i-1]!=0) {
                    for(l=0; l<n; ++l)  A[i-1-p][l] ^= A[j][l] ^= A[i-1-p][l] ^= A[j][l]; //row exchange
                    i--; //row exchanged. do it again
                    P=false;
                    break;
                }
            if(!P)  continue;
            p++;
            continue;
        }
        rank++;
        long long temp = A[i-1-p][i-1];
        A[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)  A[i-1-p][j] = (A[i-1-p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD-A[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = (A[j][k] + A[i - 1 - p][k] * mul) % MOD;
        }
    }
    return rank;
}
inline vector<vector<long long>> matrix_inverse(vector<vector<long long>> A) {
    if (A.size() != A.front().size()) {
        printf("Matrix Inversion Error : Matrix is not square\n\n");
        exit(1);
    }
    auto n = A.size();
    long long i, j, k, l;
    vector<vector<long long>> I(n, vector<long long>(n, 0));
    for (i = 0; i < n; ++i)  I[i][i] = 1;
    for (i = 1; i < n; ++i) {
        if (A[i - 1][i - 1]==0)
        {
            bool P=true;
            for(j=i; j<n; ++j)
                if(A[j][i-1]!=0)
                {
                    for(l=0; l<n; ++l)  {
                        A[j][l] ^= A[i-1][l] ^= A[j][l] ^= A[i-1][l]; //row exchange
                        I[j][l] ^= I[i-1][l] ^= I[j][l] ^= I[i-1][l];
                    }
                    P=false;
                    break;
                }
            if(P) {
                printf("Matrix Inversion Error : Matrix is singular\n\n");
                exit(1);
            }
        }
        for (j = i; j < n; ++j) {
            long long mul = (MOD-A[j][i-1]) * inverse(A[i - 1][i - 1]) % MOD;
            for (k = 0; k < n; ++k)
            {
                A[j][k] = (A[j][k] + A[i - 1][k] * mul) % MOD;
                I[j][k] = (I[j][k] + I[i - 1][k] * mul) % MOD;
            }
        }
    }
    for (i = n - 2; i >= 0; --i) {
        if (A[i + 1][i + 1]==0) {
            printf("Matrix Inversion Error : Matrix is singlular\n\n");
            exit(1);
        }
        for (j = i; j >= 0; --j) {
            long long mul = (MOD-A[j][i+1]) * inverse(A[i + 1][i + 1]) % MOD;
            for (k = 0; k < n; ++k)
            {
                A[j][k] = (A[j][k] + A[i + 1][k] * mul) % MOD;
                I[j][k] = (I[j][k] + I[i + 1][k] * mul) % MOD;
            }
        }
    }
    for (i = 0; i < n; ++i) {
        long long t = inverse(A[i][i]);
        for (j = 0; j < n; ++j)
            I[i][j] = (I[i][j] * t) % MOD;
    }
    return I;
}
inline long long matrix_determinant(vector<vector<long long>> A) {
    if (A.size() != A.front().size()) {
        printf("Matrix determinant Error : Matrix is not square\n\n");
        exit(1);
    }
    long long tr = 1;
    auto n = A.size();
    long long i, j, k, l;
    for (i = 1; i < n; ++i) {
        if (A[i - 1][i - 1]==0)
        {
            bool P=true;
            for(j=i; j<n; ++j)
                if(A[j][i-1]!=0)
                {
                    for(l=0; l<n; ++l)  A[j][l] ^= A[i-1][l] ^= A[j][l] ^= A[i-1][l]; //row exchange
                    tr*=-1;
                    P=false;
                    break;
                }
            if(P)   return 0;
        }
        for (j = i; j < n; ++j) {
            long long mul = (MOD-A[j][i - 1]) * inverse(A[i - 1][i - 1]) % MOD;
            if(!mul)  continue;
            for (k = 0; k < n; ++k)
                A[j][k] = (A[j][k] + A[i - 1][k] * mul) % MOD;
        }
    }
    long long r = tr==1?1:MOD-1;
    for(i=0; i<n; ++i)
        r = r * A[i][i] % MOD;
    return r;
}
inline vector<vector<long long>> Null_Space(vector<vector<long long>> A, bool Orth) {
    int m = (int)A.size(), n = (int)A[0].size();
    vector<int> piv;
    int i,j,k,l,p=0,rank=0;
    for(i=1; i-1<n && i-1-p<m; ++i)
    {
        if(A[i-1-p][i-1]==0) //pivot is zero
        {
            bool P=true;
            for(j=i-p; j<m; ++j) //row exchange is allowed
                if(A[j][i-1]!=0) {
                    for(l=0; l<n; ++l)  A[i-1-p][l] ^= A[j][l] ^= A[i-1-p][l] ^= A[j][l]; //row exchange
                    i--; //row exchanged. do it again
                    P=false;
                    break;
                }
            if(!P)  continue;
            p++;
            continue;
        }
        piv.push_back(i-1); //pivot location tracker
        rank++;
        long long temp = A[i-1-p][i-1];
        A[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)  A[i-1-p][j] = (A[i-1-p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD-A[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = (A[j][k] + A[i - 1 - p][k] * mul) % MOD;
        }
    }
    if(rank==n) {
        vector<vector<long long>> NR;
        return NR;
    }
    for(i=(int)piv.size()-1; i>0; --i) //upper elimination
        for(j=i-1; j>=0; --j)
        {
            long long mul = MOD-A[j][piv[i]];
            for(k=piv[i]; k<n; ++k)
                A[j][k] = (A[j][k] + A[i][k] * mul) % MOD;
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
    vector<vector<long long>> TR = matrix_transpose(A), F(n-rank, vector<long long>(rank,0));
    vector<pair<int,int>> exc;
    for(i=0; i<rank; ++i) {
        if(piv[i]!=i) {
            for(j=0; j<TR[i].size(); ++j)
                TR[i][j] ^= TR[piv[i]][j] ^= TR[i][j] ^= TR[piv[i]][j];
            exc.push_back({i,piv[i]});
            piv[i]=i;
        }
    }
    for(p=i; i<n; ++i)
        for(j=0; j<rank; ++j)
            F[i-p][j] = TR[i][j] ? MOD - TR[i][j] : 0; //remaining col of A is not from I
    vector<vector<long long>> N = matrix_transpose(F);
    for(i=0; i<F.size(); ++i) {
        vector<long long> te(F.size(),0);
        te[i]=1;
        N.push_back(te); // I padding
    }
    for(i=(int)exc.size()-1; i>=0; --i)
        for(j=0; j<N.front().size(); ++j)
            N[exc[i].first][j] ^= N[exc[i].second][j] ^= N[exc[i].first][j] ^= N[exc[i].second][j];
    if(!Orth)   return N;
    vector<vector<long long>> X = matrix_transpose(N);  //G-S process
    vector<vector<long long>> V = X;
    vector<long long> DP(2,0);
    DP[1] = (V[0]*V[0]) % MOD;
    for(i=1; i<N[0].size(); ++i)
    {
        for(j=1; j<=i; ++j)
        {
            long long c = (V[j-1] * X[i] % MOD) * inverse(DP[j]) % MOD;
            for(l=0; l<V[i].size(); ++l)
                V[i][l] = (V[i][l] + MOD - (c * V[j-1][l] % MOD)) % MOD;
        }
        DP.push_back((V[i]*V[i]) % MOD);
        if(DP.back()==0) {
            printf("NullSpace's G-S Process Error : Matrix has dependent column\n\n");
            exit(1);
        }
    }
    return matrix_transpose(V);
}
inline vector<long long> Ax_b(vector<vector<long long>>& A, vector<long long> b) {
    if(A.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }
    int m = (int)A.size(), n = (int)A[0].size();
    n++;
    vector<vector<long long>> R(m,vector<long long>(n,0));
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
                    vector<long long> temp = R[i-1-p];
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
        long long temp = R[i-1-p][i-1];
        R[i-1-p][i-1]=1;
        for(j=i; j<n; ++j)    R[i-1-p][j] = (R[i-1-p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD-R[j][i - 1];
            if(mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                R[j][k] = (R[j][k] + R[i - 1 - p][k] * mul) % MOD;
        }
    }
    for(i=(int)piv.size()-1; i>0; --i) //upper elimination
        for(j=i-1; j>=0; --j)
        {
            long long mul = MOD-R[j][piv[i]];
            for(k=piv[i]; k<n; ++k)
                R[j][k] = (R[j][k] + R[i][k] * mul) % MOD;
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
            //exit(1);
            return {0};
        }
        R.pop_back();
        b.pop_back();
    }
    if(R.size() == R[0].size()-1) {
        vector<long long> r(R.size());
        for(i=0; i<r.size(); ++i)   r[i]=R[i][n-1];
        return r;
    }
    vector<long long> r(n-1,0);
    for(i=0, p=0; i<n-1; ++i) {
        if(piv[p]==i) {
            r[i]=R[p][n-1];
            p++;
        }
    }
    return r;
}
inline void matrix_diagonalize_2x2(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    S.resize(2,vector<long long>(2,0));
    D.resize(2,vector<long long>(2,0));
    long long inroot = ((A[0][0]+A[1][1])*(A[0][0]+A[1][1])%MOD)*inverse(4)%MOD, seed=0,seed2=0;
    inroot = (inroot + (A[0][1]*A[1][0]) - ((A[0][0]*A[1][1])%MOD) + MOD)%MOD;
    if(!inroot)
        D[0][0]=D[1][1]=(A[0][0]+A[1][1])*inverse(2)%MOD;
    else {
        seed = seeds[inroot] * inverse(2) % MOD;
        seed2 = power(primitive,seed);
        long long fr = (A[0][0]+A[1][1])*inverse(2)%MOD;
        D[0][0] = (fr + seed2)%MOD;
        D[1][1] = (fr - seed2 + MOD)%MOD;
    }
    if(D[0][0]==D[1][1]) {
        S = Null_Space(A - I_n(2,D[0][0]), Orth);
        if(S.empty())
            exit(1);
        else if(S[0].size()!=2)
            exit(1);
    }
    else {
        vector<vector<long long>> ZN = Null_Space(A - I_n(2, D[0][0]), false);
        if(ZN.empty())
            exit(1);
        else if(ZN[0].size()!=1)
            exit(1);
        S[0][0]=ZN[0][0];   S[1][0]=ZN[1][0];
        ZN = Null_Space(A - I_n(2, D[1][1]), false);
        if(ZN.empty())
            exit(1);
        else if(ZN[0].size()!=1)
            exit(1);
        S[0][1]=ZN[0][0];   S[1][1]=ZN[1][0];
    }
}
inline void matrix_diagonalize_BF(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    if (A.size() != A.front().size()) {
        printf("Matrix diagonalization Error : Matrix is not square\n\n");
        exit(1);
    }
    if (!matrix_determinant(A)) {
        printf("Matrix diagonalization Error : Try Invertible matrix\n\n");
        exit(1);
    }
    if (matrix_power(A, MOD-1) != I_n((int)A.size())) {
        printf("Matrix diagonalization Error : Matrix is not diagonalizable\n\n"); //if not periodic, not diagonalizable
        exit(1);
    }
    int i,j,k,l,n=(int)A.size(), vc=0;
    S.clear();  D.clear();
    S.resize(n,vector<long long>(n,0));    D.resize(n,vector<long long>(n,0));
    vector<vector<long long>> ZN;
    long long trace = 0;
    for(i=0; i<n; ++i)  trace = (trace + A[i][i]) % MOD;
    for(i=0; i<MOD; ++i) { //eigenvalue zero to MOD-1
        ZN = Null_Space(A, Orth);
        if(!ZN.empty()) { // det(A - eigenvalue*I) == 0    <-- eigenvalue found
            for(k=0; k<ZN[0].size(); ++k)   D[k+vc][k+vc]=i;
            for(k=0; k<ZN.size(); ++k)
                for(l=0; l<ZN[0].size(); ++l)
                    S[k][l+vc]=ZN[k][l];
            vc+=(int)ZN[0].size();
            //printf(" -- total %d eigenvalues found. --> %d\n",vc,i);
            if(vc==n-1) { // only one more to go
                long long EigSum = 0;
                for(k=0; k<n-1; ++k)    EigSum = (EigSum + D[k][k]) % MOD;
                D[n-1][n-1] = (trace - EigSum + MOD) % MOD;
                //printf(" -- last eigenvalue found. --> %lld\n\n\n",D[n-1][n-1]);
                for(k=0; k<n; ++k)  A[k][k] = (A[k][k] + MOD - D[n-1][n-1] + i) % MOD;
                ZN = Null_Space(A, Orth);
                for(k=0; k<ZN.size(); ++k)  S[k][vc]=ZN[k][0];
                return;
            }
            if(vc==n)   return; //maximum n eigenvalues or eigenvectors.
        }
        for(j=0; j<n; ++j)  A[j][j] = (A[j][j] + MOD - 1) % MOD; //minus I     (we are trying every possible number as eigenvalue)
    }
    //printf("\n\n");
    //matrix_print(S);
    return;
}
inline void matrix_diagonalize_fast(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    if (A.size() != A.front().size()) {
        printf("Matrix diagonalization Error : Matrix is not square\n\n");
        exit(1);
    }
    if (!matrix_determinant(A)) {
        printf("Invertible matrix only for now\n\n");
        exit(1);
    }
    int i,j,k,l,n=(int)A.size(), vc=0, powC=(int)MOD-1, vct;
    vector<vector<long long>> ZN;
    S.resize(n,vector<long long>(n,0));
    D.resize(n,vector<long long>(n,0));
    vector<vector<long long>> PM;
    if(n==1) {
        S[0][0]=1;  D[0][0]=A[0][0];
        return;
    }
    if (matrix_power(A, MOD-1) != I_n(n)) {
        printf("Matrix diagonalization Error : Matrix is not diagonalizable\n\n"); //if not periodic, not diagonalizable
        exit(1);
    }
    if(n==2) {
        matrix_diagonalize_2x2(A, S, D, Orth);
        return;
    }
//    vector<long long> prev_roots;
//    vector<long long> roots={1,MOD-1};
//    //vector<int> eigen_count;
//    for(i=1; i<MOD_decompose.size(); ++i,vc=0) {
//        prev_roots.clear(); eigen_count.clear();
//        powC/=MOD_decompose[i-1];
//        PM = matrix_power(A, powC);
//        for(j=0; j<roots.size() && vc<n; ++j) {
//            long long rank = matrix_rank( PM - roots[j]*I );
//            if(rank!=n) {
//                vc+=n-rank;
//                prev_roots.push_back(roots[j]); //filter
//                //eigen_count.push_back((int)n-(int)rank);
//            }
//        }
//        roots.clear();
//        for(j=0; j<prev_roots.size(); ++j) {     //most time consuming part (was)
//            long long seed = seeds[prev_roots[j]] * inverse(MOD_decompose[i]) % MOD;
//            long long seed2 = power(primitive,seed);
//            for(l=0; l<ones_roots[MOD_decompose[i]].size(); ++l)
//                roots.push_back(seed2 * ones_roots[MOD_decompose[i]][l] % MOD);  // in MOD_decompose, there is no p-1. max is (p-1)/2
//        }
//    }
    
    vector<long long> prev_roots;
    vector<vector<long long>> roots(1,vector<long long>());
    vector<long long> eigen_count(1,n), ect;
    for(i=0; i<ones_roots[MOD_decompose.back()].size(); ++i)
        roots[0].push_back(ones_roots[MOD_decompose.back()][i]);
    for(i=(int)MOD_decompose.size()-2; i>=0; --i,vc=0,prev_roots.clear(),ect.clear()) {
        powC/=MOD_decompose[i+1];
        PM = matrix_power(A, powC);
        for(k=0; k<roots.size(); ++k) {
            for(j=vc=0; j<roots[k].size() && vc<eigen_count[k]; ++j) {
                long long rank = matrix_rank( PM - I_n(n,roots[k][j]) ); // test
                if(rank!=n) {
                    vc+=n-rank;
                    prev_roots.push_back(roots[k][j]);
                    ect.push_back(n-rank);
                }
            }
        }
        roots.clear();
        roots.resize(prev_roots.size(), vector<long long>());
        eigen_count=ect;
        for(j=0; j<prev_roots.size(); ++j) {     //most time consuming part (was)
            long long seed = seeds[prev_roots[j]] * inverse(MOD_decompose[i]) % MOD;
            long long seed2 = power(primitive,seed);
            for(l=0; l<ones_roots[MOD_decompose[i]].size(); ++l)
                roots[j].push_back(seed2 * ones_roots[MOD_decompose[i]][l] % MOD);  // in MOD_decompose, there is no p-1. max is (p-1)/2
        }
    }
    
    long long trace = 0;
    for(i=0; i<n; ++i)  trace = (trace + A[i][i]) % MOD;
    for(k=vc=0; k<roots.size(); ++k) {
        for(j=vct=0; j<2 && vct<eigen_count[k]; ++j) {   //roots[?].size() == 2
            ZN = Null_Space(A - I_n(n,roots[k][j]), Orth);
            if(!ZN.empty()) {
                for(i=0; i<ZN[0].size(); ++i) {
                    D[i+vc][i+vc]=roots[k][j];
                    for(l=0; l<ZN.size(); ++l)
                        S[l][i+vc]=ZN[l][i];
                }
                vc+=(int)ZN[0].size();
                vct+=(int)ZN[0].size();
                if(vc==n-1) {
                    long long EigSum = 0;
                    for(i=0; i<n-1; ++i)    EigSum = (EigSum + D[i][i]) % MOD;
                    D[n-1][n-1] = (trace - EigSum + MOD) % MOD;
                    //printf(" -- last eigenvalue found. --> %lld\n\n\n",D[n-1][n-1]);
                    ZN = Null_Space(A - I_n(n,D[n-1][n-1]), Orth);
                    for(i=0; i<ZN.size(); ++i)  S[i][vc]=ZN[i][0];
                    return;
                }
            }
        }
    }
    return;
}
inline void matrix_diagonalize_fast2(vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth, long long powC, int pi, long long ev) {
    int i,j,k,l,n=(int)A.size(), vc=0, vct;
    vector<vector<long long>> ZN;
    S.resize(n,vector<long long>(n,0));
    D.resize(n,vector<long long>(n,0));
    vector<vector<long long>> PM;
    if(n==2) {
        matrix_diagonalize_2x2(A, S, D, Orth);
        return;
    }
    
    vector<long long> prev_roots; //(ev)
    vector<vector<long long>> roots(1,vector<long long>());
    vector<long long> eigen_count(1,n), ect;
    long long seed = seeds[ev] * inverse(MOD_decompose[pi]) % MOD;
    long long seed2 = power(primitive,seed);
    for(i=0; i<ones_roots[MOD_decompose[pi]].size(); ++i)
        roots[0].push_back(seed2 * ones_roots[MOD_decompose[pi]][i] % MOD);
    for(i=pi+1; i<MOD_decompose.size(); ++i,vc=0,prev_roots.clear(),ect.clear()) {
        powC/=MOD_decompose[i-1];
        
        PM = matrix_power(A, powC);
        for(k=0; k<roots.size(); ++k) {
            for(j=vc=0; j<roots[k].size() && vc<eigen_count[k]; ++j) {
                long long rank = matrix_rank( PM - I_n(n,roots[k][j]) ); // test
                if(rank!=n) {
                    vc+=n-rank;
                    prev_roots.push_back(roots[k][j]);
                    ect.push_back(n-rank);
                }
            }
        }
        
        roots.clear();
        roots.resize(prev_roots.size(), vector<long long>());
        eigen_count=ect;
        for(j=0; j<prev_roots.size(); ++j) {     //prev_roots.size() is max n
            seed = seeds[prev_roots[j]] * inverse(MOD_decompose[i]) % MOD;
            seed2 = power(primitive,seed);
            for(l=0; l<ones_roots[MOD_decompose[i]].size(); ++l)
                roots[j].push_back(seed2 * ones_roots[MOD_decompose[i]][l] % MOD);  // in MOD_decompose, there is no p-1. max is (p-1)/2
        }
    }
    
    long long trace = 0;
    for(i=0; i<n; ++i)  trace = (trace + A[i][i]) % MOD;
    for(k=vc=0; k<roots.size(); ++k) {
        for(j=vct=0; j<roots[k].size() && vct<eigen_count[k]; ++j) {
            ZN = Null_Space(A - I_n(n,roots[k][j]), Orth);
            if(!ZN.empty()) {
                for(i=0; i<ZN[0].size(); ++i) {
                    D[i+vc][i+vc]=roots[k][j];
                    for(l=0; l<ZN.size(); ++l)
                        S[l][i+vc]=ZN[l][i];
                }
                vc+=(int)ZN[0].size();
                vct+=(int)ZN[0].size();
                if(vc==n-1) {
                    long long EigSum = 0;
                    for(i=0; i<n-1; ++i)    EigSum = (EigSum + D[i][i]) % MOD;
                    D[n-1][n-1] = (trace - EigSum + MOD) % MOD;
                    //printf(" -- last eigenvalue found. --> %lld\n\n\n",D[n-1][n-1]);
                    ZN = Null_Space(A - I_n(n,D[n-1][n-1]), Orth);
                    for(i=0; i<ZN.size(); ++i)  S[i][vc]=ZN[i][0];
                    return;
                }
            }
        }
    }
    return;
}
inline void matrix_diagonalize_henry(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    if (A.size() != A.front().size()) {
        printf("Matrix diagonalization Error : Matrix is not square\n\n");
        exit(1);
    }
    if (!matrix_determinant(A)) {
        printf("Invertible matrix only for now\n\n");
        exit(1);
    }
    int i,j,k,n=(int)A.size(), vc=0, mati=0;
    vector<vector<long long>> ZN;
    S = I_n(n);
    D.clear();
    if (matrix_power(A, MOD-1) != I_n(n)) {
        printf("Matrix diagonalization Error : Matrix is not diagonalizable\n\n"); //if not periodic, not diagonalizable
        exit(1);
    }
    vector<vector<vector<long long>>> M;    M.push_back(A);
    vector<long long> FE(1,1);
    long long powC=MOD-1;
    int pi=0;
    for(int stp=0; pi<MOD_decompose.size()>>1; ++pi,stp=0) {
    //for(int pi=0,stp=0; pi<5; ++pi,stp=0) {
        int matiu = (int)M.size();
        powC/=MOD_decompose[pi];
        vector<vector<long long>> ST(n, vector<long long>(n,0));
        for(; mati<matiu; ++mati,vc=0) {
            if(M[mati].size()==1) {
                M.push_back(M[mati]);
                FE.push_back(M[mati][0][0]);
                ST[stp][stp]=1;
                stp++;
                continue;
            }
            vector<vector<long long>> St(M[mati].size(), vector<long long>(M[mati].size()));
            vector<vector<long long>> PM = matrix_power(M[mati], powC);
            vector<int> esc;
            
            long long seed = seeds[FE[mati]] * inverse(MOD_decompose[pi]) % MOD;
            long long seed2 = power(primitive,seed);
            for(i=0; i<ones_roots[MOD_decompose[pi]].size(); ++i) {
                //roots[j].push_back(seed2 * ones_roots[MOD_decompose[i]][l] % MOD);  // in MOD_decompose, there is no p-1. max is (p-1)/2
            //for(i=0; i<ones_roots[(MOD-1)/powC].size() && vc<M[mati].size(); ++i) { //create S
                //ZN = Null_Space(PM - I_n((int)M[mati].size(),ones_roots[(MOD-1)/powC][i]), Orth);
                ZN = Null_Space(PM - I_n((int)M[mati].size(), seed2 * ones_roots[MOD_decompose[pi]][i] % MOD), Orth);
                if(ZN.empty())  continue;
                esc.push_back((int)ZN[0].size());
                //FE.push_back(ones_roots[(MOD-1)/powC][i]);
                FE.push_back(seed2 * ones_roots[MOD_decompose[pi]][i] % MOD);
                for(j=0; j<ZN.size(); ++j)
                    for(k=0; k<ZN[0].size(); ++k)
                        St[j][k+vc]=ZN[j][k];
                vc+=(int)ZN[0].size();
            }
            
            
            vector<vector<long long>> mt = matrix_inverse(St) * M[mati] * St;
            for(i=0; i<St.size(); ++i)
                for(j=0; j<St.size(); ++j)
                    ST[i+stp][j+stp]=St[i][j];
            stp+=St.size();
            matrix_chop(M, mt, esc);
        }
        S = S*ST;
    }
    vector<vector<long long>> ST;
    for(; mati<M.size(); ++mati) {
        if(M[mati].size()==1) {
            D = D|M[mati];
            vector<vector<long long>> St(1, vector<long long>(1,1));
            ST = ST|St;
            continue;
        }
        vector<vector<long long>> St,Dt;
        //matrix_diagonalize_fast2(M[mati], St, Dt, Orth, powC,pi,FE[mati]);
        matrix_diagonalize_fast(M[mati], St, Dt, Orth);
        ST = ST|St;
        D = D|Dt;
    }
    S = S*ST;
}

inline void func1() {
    int N=300,i,j,k;
    double avt=0;
    vector<vector<long long>> I(N,vector<long long>(N,0)),S1,D1,S2,D2;
    for(i=0; i<N; ++i)  I[i][i]=1;
    for(int trial=1; trial<1000000000; ++trial) {
        vector<vector<long long>> tm = I;
        for(i=0; i<N-1; ++i) {
            for(j=i+1; j<N; ++j) {
                long long mul = rand()%MOD;
                for(k=0; k<N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        for(i=N-1; i>0; --i) {
            for(j=i-1; j>=0; --j) {
                long long mul = rand()%MOD;
                for(k=0; k<N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        vector<vector<long long>> E(N, vector<long long>(N,0));
        for(i=0; i<N; ++i)  E[i][i] = rand()%(MOD-1)+1;
        vector<vector<long long>> DC = tm * E * matrix_inverse(tm);
        clock_t s2,f2;
        s2=clock();
        matrix_diagonalize_henry(DC, S2, D2, false);
        f2=clock();
        if(DC * S2 != S2 * D2 || matrix_determinant(S2)==0) {
            printf("NOT GOOD...\n\n");
            matrix_print(DC);
            printf("WRONG : \n");
            matrix_print(S2);
            matrix_print(D2);
            matrix_print(S2 * D2 * matrix_inverse(S2));
            exit(1);
        }
        double d2=(double)(f2-s2)/CLOCKS_PER_SEC;
        avt+=d2;
        printf("-- %d\t\t%lf sec.\t\t(avg %lf sec)\n",trial,d2,avt/trial);
    }
}
inline void func2() {
    long long po = 9223372036854775807; // 9223372036854775807;
    long long modpo = po%(MOD-1);
    vector<vector<long long>> A = {
        {63,63, 0,48,13},
        {48, 5,89,65,57},
        {32,69, 1,57,68},
        {95, 8,46,53,32},
        {34,100,50,80,70}
    },D,S;
    matrix_diagonalize_fast(A, S, D, false);
    for(int i=0; i<A.size(); ++i)  D[i][i] = power(D[i][i],modpo);
    if(matrix_power(A, po) == S * D * matrix_inverse(S))    printf("GOOD\n\n");
    else    printf("NOOOOT GOOD....\n\n");
    po=0;
}
inline void func3() {
    long long i,c1=0,cp_1=0;
    //vector<long long> a,b;
    for(i=1; i<=MOD>>1; ++i) {
        long long k = power(i,(MOD-1)>>1);
        printf("%lld\t%lld\n",i,k);
        if(k==1)    {
            //a.push_back(i);
            c1++;
        }
        else    {
            //b.push_back(i);
            cp_1++;
        }
    }
    i=0;
    return;
}
inline void func4() {
    int N=20,i,j,k;
    double avt=0;
    vector<vector<long long>> I,S1,D1,S2,D2;
    I=I_n(N);
    for(int trial=1; trial<1000000000; ++trial) {
        vector<vector<long long>> tm = I;
        for(i=0; i<N-1; ++i) {
            for(j=i+1; j<N; ++j) {
                long long mul = rand()%MOD;
                for(k=0; k<N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        for(i=N-1; i>0; --i) {
            for(j=i-1; j>=0; --j) {
                long long mul = rand()%MOD;
                for(k=0; k<N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        vector<vector<long long>> E(N, vector<long long>(N,0));
        for(i=0; i<N; ++i)  E[i][i] = rand()%(MOD-1)+1;
        vector<vector<long long>> DC = tm * E * matrix_inverse(tm);
        clock_t s1,f1,s2,f2;
        s1=clock();
        matrix_diagonalize_fast(DC, S1, D1, false);
        f1=clock();
        s2=clock();
        matrix_diagonalize_henry(DC, S2, D2, false);
        f2=clock();
        if(DC * S2 != S2 * D2 || matrix_determinant(S2)==0) {
            printf("NOT GOOD...\n\n");
            matrix_print(DC);
            printf("RIGHT : \n");
            matrix_print(S1);
            matrix_print(D1);
            matrix_print(S1 * D1 * matrix_inverse(S1));
            printf("\n\n\nWRONG : \n");
            matrix_print(S2);
            matrix_print(D2);
            matrix_print(S2 * D2 * matrix_inverse(S2));
            exit(1);
        }
        double d1=(double)(f1-s1)/CLOCKS_PER_SEC, d2=(double)(f2-s2)/CLOCKS_PER_SEC;
        avt+=d1/d2;
        printf("-- %d\t\t%lf sec vs %lf sec.\t\t%lf time faster!!\t\t(avg %lf time faster)\n",trial,d1,d2,d1/d2,avt/trial);
        //matrix_print(tm);
    }
}

int main()
{
    MOD = 100000007;
    //MOD = 131071;
    //MOD = 524287;
    //MOD = 65537;
    //MOD = 653659;
    //MOD = 101;
    Initiation();
    func4();
    vector<vector<long long>> A = {
        //        {3,5,7,2},
        //        {1,4,7,2},
        //        {6,3,9,17},
        //        {13,5,4,16}
        
        //        {1,0,0,0},
        //        {0,2,9,0},
        //        {0,9,3,1},
        //        {0,0,1,4}
        
        //        {1,2,3,4},
        //        {2,3,4,5},
        //        {9,4,5,6},
        //        {1,5,2,7}
        
        //        {3,5,9,3},
        //        {1,0,0,0},
        //        {0,1,0,0},
        //        {0,0,1,0}
        
        //        {1,2,0},
        //        {2,3,4},
        //        {3,6,1}
        
//        {63,63, 0,48,13},
//        {48, 5,89,65,57},
//        {32,69, 1,57,68},
//        {95, 8,46,53,32},
//        {34,100,50,80,70}
        
//        {{51 ,   78   , 50  ,  8  ,  42},
//            {32 ,   15   , 17   , 68 ,   47},
//            {11  ,  80 ,   77  ,  77   , 94},
//            {11   , 53  ,  36  ,  88   , 12},
//            {65 ,   77   , 8  ,  58   , 31}}
        
        {{98   , 50  ,  4 ,   80 ,   51},
            {   27  ,  68 ,   96   , 84  ,  82},
            {64  ,  57  ,  20  ,  9  ,  48},
            {     0   , 40  ,  10  ,  34   , 5},
            { 23   , 73   , 80   , 43   , 80}}
        
        
//        {1,0,0,0},
//        {0,0,0,0},
//        {0,0,2,0},
//        {0,0,0,3}
    },S,D,S1,D1,M1,M2,M3;
    vector<vector<long long>> E = {
        {1,0,0,0,0},
        {0,9,0,0,0},
        {0,0,19,0,0},
        {0,0,0,6685,0},
        {0,0,0,0,2347823}
        
//        {1,0,0},
//        {0,4,0},
//        {0,0,6}
    };
    //vector<vector<long long>> MT = A*E*matrix_inverse(A);
    
    //matrix_print(MT);
    matrix_diagonalize_fast(A, S, D, false);
    //matrix_print(D);    matrix_print(S);    matrix_print(S*D*matrix_inverse(S));  printf("\n\n");
    matrix_diagonalize_henry(A, S1, D1, false);
    matrix_print(D1);    matrix_print(S1);    matrix_print(S1*D1*matrix_inverse(S1));
    return 0;
}
