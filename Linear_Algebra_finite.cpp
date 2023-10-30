#include <vector>
#include <iostream>
using namespace std;
long long MOD = 10e9; //should be inside of INT32 range. It's GF so MOD must be a prime.
vector<int> int_inverse;
void Initiation() {
    int_inverse.resize(MOD,0);
}
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
            R[i][j] = (a[i][j] - b[i][j] + MOD) % MOD;
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
inline vector<vector<long long>> Null_Space(vector<vector<long long>> A) {
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
            for(k=piv[i]; k<n; ++k) {
                A[j][k] += (A[i][k] * mul);
                A[j][k] %= MOD;
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
        A.pop_back();
    }
    vector<vector<long long>> TR = matrix_transpose(A), F(n-rank, vector<long long>(rank,0));
    vector<pair<int,int>> exc;
    for(i=0; i<rank; ++i) {
        if(piv[i]!=i) {
            vector<long long> te = TR[i];
            TR[i] = TR[piv[i]];
            TR[piv[i]] = te;
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
    for(i=(int)exc.size()-1; i>=0; --i) {
        vector<long long> te = N[exc[i].first];
        N[exc[i].first] = N[exc[i].second];
        N[exc[i].second] = te;
    }
    return N;
}
inline void matrix_diagonalize(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D) {
    if (A.size() != A.front().size()) {
        printf("Matrix determinant Error : Matrix is not square\n\n");
        exit(1);
    }
    int i,j,k,l,n=(int)A.size(), vc=0,c=0,cp=1, dMOD = (int)MOD/100;
    S.resize(n,vector<long long>(n,0));    D.resize(n,vector<long long>(n,0));
    vector<vector<long long>> ZN;
    long long trace = 0;
    for(i=0; i<n; ++i)  trace = (trace + A[i][i]) % MOD;
    for(i=0; i<MOD; ++i,++c) { //eigenvalue zero to MOD-1
        if(c==dMOD) {
            printf(" -- %d%%.\n",cp++);
            c=0;
        }
        ZN = Null_Space(A);
        if(!ZN.empty()) { // det(A - eigenvalue*I) == 0    <-- eigenvalue found
            for(k=0; k<ZN[0].size(); ++k)   D[k+vc][k+vc]=i;
            for(k=0; k<ZN.size(); ++k)
                for(l=0; l<ZN[0].size(); ++l)
                    S[k][l+vc]=ZN[k][l];
            vc+=(int)ZN[0].size();
            printf(" -- %d eigenvalue found. --> %d\n",vc,i);
            if(vc==n-1) { // only one more to go
                long long EigSum = 0;
                for(k=0; k<n-1; ++k)    EigSum = (EigSum + D[k][k]) % MOD;
                D[n-1][n-1] = (trace - EigSum + MOD) % MOD;
                printf(" -- last eigenvalue found. --> %lld\n",D[n-1][n-1]);
                for(k=0; k<n; ++k)  A[k][k] = (A[k][k] + MOD - D[n-1][n-1] + i) % MOD;
                ZN = Null_Space(A);
                for(k=0; k<ZN.size(); ++k)  S[k][vc]=ZN[k][0];
                return;
            }
            if(vc==n)   return; //maximum n eigenvalues or eigenvectors.
        }
        for(j=0; j<n; ++j)  A[j][j] = (A[j][j] + MOD - 1) % MOD; //minus one
    }
    printf("\n\n");
}

int main()
{
//    long long i,j,t;
//    for (i = 3; i <= 100000007; i+=2) {
//        t = sqrt(i);
//        bool P=true;
//        for (j = 2; j <= t; ++j) {
//            if (!(i % j)) {
//                P=false;
//                break;
//            }
//        }
//        if(!P)  continue;
//        //
//        printf("prime : %lld\n",i);
//        vector<bool> v(i,false);
//        for (j=0; j<i; ++j) {
//            t = (j*j)%i;
//            if(v[t]) {
//                P=false;
//                break;
//            }
//            v[t]=true;
//        }
//        if(P) {
//            printf("POSSIBLE\n",i);
//            continue;
//        }
//    }
    
    
    //MOD = 100000007;
    MOD = 7;
    Initiation();
    vector<vector<long long>> A = {
//        {3,5,7,2},
//        {1,4,7,2},
//        {6,3,9,17},
//        {13,5,4,16}
        
//        {1,0,0,0},
//        {0,2,0,0},
//        {0,0,4,0},
//        {0,0,0,8}
        
        {2,0},
        {0,4}
        
//        {1,1},
//        {1,0}
        
//        {1,2,3,4},
//        {2,3,4,5},
//        {9,4,5,6},
//        {1,5,6,7}
        
//        {1,0,0,1},
//        {0,2,0,0},
//        {0,0,3,0},
//        {1,1,0,4}
        
//        {1,0,0,0},
//        {0,0,0,0},
//        {0,0,1,0},
//        {0,0,0,1}
    },S,D;
    
//    vector<vector<long long>> NS = Null_Space(A);
//    matrix_print(NS);
//    matrix_print(A * NS);
    
    matrix_diagonalize(A, S, D);
    matrix_print(S);
    matrix_print(D);
    vector<vector<long long>> A2 = S * D * matrix_inverse(S);
    matrix_print(A2);
    
    //long long po = 9223372036854775807;
    long long po = 3;
    for(int i=0; i<A.size(); ++i)  D[i][i] = power(D[i][i],po);
    matrix_print(S*D*matrix_inverse(S));
    if(matrix_power(A, po) == S * D * matrix_inverse(S))
        printf("GOOD\n\n");
    
    return 0;
}