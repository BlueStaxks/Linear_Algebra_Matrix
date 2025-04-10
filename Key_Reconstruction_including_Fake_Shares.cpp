#include <fstream>
#include <time.h>
#include <chrono>

#include <set>
#include <omp.h>
#include <math.h>
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
    if (!a) {
        printf("Integer Inverse Error : 0 has no inverse.\n\n");
        exit(1);
    }
    if (int_inverse[a])  return int_inverse[a];
    long long q, r1 = MOD, r2 = a, r = 1, t1 = 0, t2 = 1, t;
    while (r) {
        q = r1 / r2;    r = r1 % r2;
        t = (t1 - q * t2) % MOD;
        if (t < 0) t += MOD;
        r1 = r2;  r2 = r;   t1 = t2;  t2 = t;
    }
    //    if(r1!=1) {
    //        printf("Integer Inverse Error : %lld does not have inverse in Z(%lld)\n\n",a,MOD);
    //        exit(1);
    //    }
    int_inverse[a] = (int)t1;
    return t1;
}
inline long long power(long long a, long long n) {
    long long y = 1;
    while (n > 0) {
        if (n & 1)
            y = (a * y) % MOD;
        a = (a * a) % MOD;
        n >>= 1;
    }
    return y;
}
inline long long gcd(long long a, long long b) {
    if (a < b) a ^= b ^= a ^= b;
    long long n;
    while (b) {
        n = a % b;
        a = b;
        b = n;
    }
    return a;
}
inline vector<long long> decompose(long long a) {
    vector<long long> r;
    while (!(a & 1)) {
        r.push_back(2);
        a >>= 1;
    }
    for (long long i = 3; i <= sqrt(a); i += 2)
        if (a % i == 0) {
            r.push_back(i);
            a /= i;
        }
    r.push_back(a);
    return r;
}
inline vector<long long> divisor(long long a) {
    vector<long long> r;
    long long sq = sqrt(a);
    for (long long i = 1; i <= sq; ++i)
        if (a % i == 0) {
            r.push_back(i);
            if (i * i != a) r.push_back(a / i);
        }
    sort(r.begin(), r.end());
    return r;
}
void readData(const string& filename) {
    // Open file for binary reading
    ifstream inFile(filename, ios::binary);
    if (!inFile) {
        cerr << "Error opening file for reading." << endl;
        return;
    }
    // Function to read vectors from file
    auto readVector = [&inFile](auto& vec) {
        size_t size;
        inFile.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(decltype(vec[0])));
    };
    // Read the long long primitive
    inFile.read(reinterpret_cast<char*>(&primitive), sizeof(primitive));
    // Read vectors
    readVector(int_inverse);
    readVector(seeds);
    readVector(MOD_decompose);
    readVector(MOD_divisors);
    // Read vector of vectors<int>
    size_t outerSize;
    inFile.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));
    ones_roots.resize(outerSize);
    for (auto& vec : ones_roots) {
        readVector(vec);
    }
    inFile.close();
    // Check for errors
    if (!inFile.good()) {
        cerr << "Error occurred during file read." << endl;
        return;
    }
}
void writeData(const string& filename) {
    // Open a file for binary writing, ios::binary tells ofstream to treat the file as binary
    ofstream outFile(filename, ios::binary | ios::out);
    if (!outFile) {
        cerr << "Error opening file for writing." << endl;
        return;
    }
    // Write the long long primitive
    outFile.write(reinterpret_cast<char*>(&primitive), sizeof(primitive));
    // Lambda to write vector to file, handles both vector<int> and vector<long long>
    auto writeVector = [&outFile](const auto& vec) {
        size_t size = vec.size();
        outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
        outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(decltype(vec[0])));
    };
    // Write vector<int> and vector<long long> to file
    writeVector(int_inverse);
    writeVector(seeds);
    writeVector(MOD_decompose);
    writeVector(MOD_divisors);
    // Write vector of vectors<int>
    size_t outerSize = ones_roots.size();
    outFile.write(reinterpret_cast<const char*>(&outerSize), sizeof(outerSize));
    for (const auto& vec : ones_roots) {
        writeVector(vec);
    }
    outFile.close();
    if (!outFile.good()) {
        cerr << "Error occurred during file write." << endl;
        return;
    }
    //cout << "Data written successfully." << endl;
}
void Initiation() {
    string name = to_string(MOD);
    name+=".bin";
    ifstream file1(name);
    if(!file1) { //no file
        int_inverse.resize(MOD, 0);
        MOD_decompose = decompose(MOD - 1);
        MOD_divisors = divisor(MOD - 1);
        primitive = 0;
        for (int i = 2; i < MOD; ++i) {
            bool P = true;
            for (int j = 1; j < MOD_divisors.size() - 1; ++j)
                if (power(i, MOD_divisors[j]) == 1) {
                    P = false;
                    break;
                }
            if (P) {
                primitive = i; //find smallest primitive root
                break;
            }
        }
        seeds.resize(MOD);
        ones_roots.resize(MOD);
        for (long long i = 1, t = primitive; i < MOD; ++i, t = t * primitive % MOD)
            seeds[t] = (int)i;  //primitive ^ seeds[i] = i
        #pragma omp parallel for schedule(dynamic)
        for (long long i = 1; i < MOD; ++i) {
            for (int j = 0; j < MOD_divisors.size() - 1; ++j) // 1^(1/i) = ones_roots.front() ~ back()
                if (power(i, MOD_divisors[j]) == 1) {
                    #pragma omp critical
                    ones_roots[MOD_divisors[j]].push_back((int)i); //calculate order of all numbers
                }
        } //files are created.
        writeData(name);
    }
    else
        readData(name);
}

template <typename T>
inline bool operator == (vector<vector<T>>& a, vector<vector<T>>& b) {
    if (a.front().size() != b.front().size() || a.size() != b.size())  return false;
    for (int i = 0; i < a.size(); ++i)
        for (int j = 0; j < a[0].size(); ++j)
            if (a[i][j] != b[i][j])
                return false;
    return true;
}
template <typename T>
inline bool operator != (vector<vector<T>>& a, vector<vector<T>>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator == (vector<T>& a, vector<T>& b) {
    if (a.size() != b.size())  return false;
    for (int i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}
template <typename T>
inline bool operator != (vector<T>& a, vector<T>& b) {
    return !(a == b);
}
inline void matrix_print(vector<vector<long long>> a) {
    if (a.empty())   return;   if (a[0].empty())    return;
    for (int i = 0; i < a.size(); ++i)
    {
        for (int j = 0; j < a[0].size(); ++j)
            printf("%lld\t", a[i][j]);
        printf("\n");
    }
    printf("\n\n");
}
inline void vector_print(vector<long long> a) {
    for (int i = 0; i < a.size(); ++i)
        printf("%lld\t", a[i]);
    printf("\n\n\n");
}


inline vector<vector<long long>> operator * (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.front().size() != b.size()) {
        printf("Matrix Multiplication Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(), 0));
    int i, j, k;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            for (k = 0; k < b.size(); ++k)
                R[i][j] = (R[i][j] + a[i][k] * b[k][j]) % MOD;
    return R;
}
inline vector<long long> operator * (const vector<vector<long long>>& a, const vector<long long>& b) {
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
inline long long operator * (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Dot Product Error : Vector size does not match\n\n");
        exit(1);
    }
    long long r = 0;
    for (auto i = 0; i < a.size(); ++i)
        r = (r + a[i] * b[i]) % MOD;
    return r;
}
inline vector<long long> operator * (const long long& a, const vector<long long>& b) {
    vector<long long> R(b.size());
    for (auto i = 0; i < b.size(); ++i)
        R[i] = (a * b[i]) % MOD;;
    return R;
}
inline vector<vector<long long>> operator * (const long long& a, const vector<vector<long long>>& b) {
    vector<vector<long long>> R(b.size(), vector<long long>(b.front().size()));
    for (auto i = 0; i < b.size(); ++i)
        for (auto j = 0; j < b.front().size(); ++j)
            R[i][j] = (a * b[i][j]) % MOD;;
    return R;
}
inline vector<vector<long long>> operator + (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Addition Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(), 0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = (a[i][j] + b[i][j]) % MOD;
    return R;
}
inline vector<vector<long long>> operator - (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.front().size() != b.front().size() || a.size() != b.size()) {
        printf("Matrix Subtraction Error : Matrix size does not match\n\n");
        exit(1);
    }
    vector<vector<long long>> R(a.size(), vector<long long>(b.front().size(), 0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < b.front().size(); ++j)
            R[i][j] = (a[i][j] - b[i][j] + MOD) % MOD;
    return R;
}
inline vector<vector<long long>> operator | (const vector<vector<long long>>& a, const vector<vector<long long>>& b) { //diagonal expansion
    if (a.empty())   return b;
    if (b.empty())   return a;
    vector<vector<long long>> R(a.size() + b.size(), vector<long long>(a.front().size() + b.front().size(), 0));
    int i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < a.front().size(); ++j)
            R[i][j] = a[i][j];
    for (i = (int)a.size(); i < R.size(); ++i)
        for (j = (int)a.front().size(); j < R.front().size(); ++j)
            R[i][j] = b[i - a.size()][j - a.front().size()];
    return R;
}
inline vector<long long> operator + (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Addition Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<long long> R(a.size(), 0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = (a[i] + b[i]) % MOD;
    return R;
}
inline vector<long long> operator - (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Subtraction Error : Vector size does not match\n\n");
        exit(1);
    }
    vector<long long> R(a.size(), 0);
    for (int i = 0; i < a.size(); ++i)
        R[i] = (a[i] - b[i]) % MOD;
    return R;
}
inline vector<long long> Extended_Euclid(long long a, long long b) {
    long long q, r1 = a, r2 = b, r = 1, s1 = 1, s2 = 0, s, t1 = 0, t2 = 1, t;
    while (r) {
        q = r1 / r2;    r = r1 % r2;
        s = (s1 - q * s2);
        t = (t1 - q * t2);
        r1 = r2;  r2 = r;   s1 = s2;  s2 = s;   t1 = t2;  t2 = t;
    }
    return { r1,s1,t1 };
}
inline vector<vector<long long>> I_n(int n) {
    vector<vector<long long>> I(n, vector<long long>(n, 0));
    for (int i = 0; i < n; ++i)  I[i][i] = 1;
    return I;
}
inline vector<vector<long long>> I_n(int n, long long a) {
    vector<vector<long long>> I(n, vector<long long>(n, 0));
    for (int i = 0; i < n; ++i)  I[i][i] = a;
    return I;
}
inline vector<vector<long long>> matrix_transpose(const vector<vector<long long>>& a) {
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
    for (int i = 0, p = 0; i < list.size(); ++i) {
        M.push_back(N);
        M.back().resize(list[i], vector<long long>(list[i]));
        for (int j = 0; j < list[i]; ++j)
            for (int k = 0; k < list[i]; ++k)
                M.back()[j][k] = F[j + p][k + p];
        p += list[i];
    }
}
inline vector<vector<long long>> matrix_partial_multiply(vector<vector<long long>>& A, vector<vector<long long>>& B, vector<int> list) {
    vector<vector<long long>> R(A.size(), vector<long long>(A.size(), 0));
    int i, j, l, k, n = (int)A.size(), p = 0;
    for (i = 0; i < list.size(); ++i) {
        for (j = 0; j < n; ++j)
            for (l = p; l < p + list[i]; ++l)
                for (k = 0; k < list[i]; ++k)
                    R[j][l] = (R[j][l] + A[j][p + k] * B[p + k][l]) % MOD;
        p += list[i];
    }
    return R;
}
inline long long matrix_rank(vector<vector<long long>> A) {
    int m = (int)A.size(), n = (int)A[0].size();
    int i, j, k, l, p = 0;
    long long rank = 0;
    for (i = 1; i - 1 < n && i - 1 - p < m; ++i)
    {
        if (A[i - 1 - p][i - 1] == 0) //pivot is zero
        {
            bool P = true;
            for (j = i - p; j < m; ++j) //row exchange is allowed
                if (A[j][i - 1] != 0) {
                    for (l = 0; l < n; ++l)  A[i - 1 - p][l] ^= A[j][l] ^= A[i - 1 - p][l] ^= A[j][l]; //row exchange
                    i--; //row exchanged. do it again
                    P = false;
                    break;
                }
            if (!P)  continue;
            p++;
            continue;
        }
        rank++;
        long long temp = A[i - 1 - p][i - 1];
        A[i - 1 - p][i - 1] = 1;
        for (j = i; j < n; ++j)  A[i - 1 - p][j] = (A[i - 1 - p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD - A[j][i - 1];
            if (mul == 0)  continue;
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
        if (A[i - 1][i - 1] == 0)
        {
            bool P = true;
            for (j = i; j < n; ++j)
                if (A[j][i - 1] != 0)
                {
                    for (l = 0; l < n; ++l) {
                        A[j][l] ^= A[i - 1][l] ^= A[j][l] ^= A[i - 1][l]; //row exchange
                        I[j][l] ^= I[i - 1][l] ^= I[j][l] ^= I[i - 1][l];
                    }
                    P = false;
                    break;
                }
            if (P) {
                printf("Matrix Inversion Error : Matrix is singular\n\n");
                exit(1);
            }
        }
        for (j = i; j < n; ++j) {
            long long mul = (MOD - A[j][i - 1]) * inverse(A[i - 1][i - 1]) % MOD;
            for (k = 0; k < n; ++k)
            {
                A[j][k] = (A[j][k] + A[i - 1][k] * mul) % MOD;
                I[j][k] = (I[j][k] + I[i - 1][k] * mul) % MOD;
            }
        }
    }
    for (i = n - 2; i >= 0; --i) {
        if (A[i + 1][i + 1] == 0) {
            printf("Matrix Inversion Error : Matrix is singlular\n\n");
            exit(1);
        }
        for (j = i; j >= 0; --j) {
            long long mul = (MOD - A[j][i + 1]) * inverse(A[i + 1][i + 1]) % MOD;
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
        if (A[i - 1][i - 1] == 0)
        {
            bool P = true;
            for (j = i; j < n; ++j)
                if (A[j][i - 1] != 0)
                {
                    for (l = 0; l < n; ++l)  A[j][l] ^= A[i - 1][l] ^= A[j][l] ^= A[i - 1][l]; //row exchange
                    tr *= -1;
                    P = false;
                    break;
                }
            if (P)   return 0;
        }
        for (j = i; j < n; ++j) {
            long long mul = (MOD - A[j][i - 1]) * inverse(A[i - 1][i - 1]) % MOD;
            if (!mul)  continue;
            for (k = 0; k < n; ++k)
                A[j][k] = (A[j][k] + A[i - 1][k] * mul) % MOD;
        }
    }
    long long r = tr == 1 ? 1 : MOD - 1;
    for (i = 0; i < n; ++i)
        r = r * A[i][i] % MOD;
    return r;
}
inline vector<vector<long long>> Null_Space(vector<vector<long long>> A, bool Orth) {
    int m = (int)A.size(), n = (int)A[0].size();
    vector<int> piv;
    int i, j, k, l, p = 0, rank = 0;
    for (i = 1; i - 1 < n && i - 1 - p < m; ++i)
    {
        if (A[i - 1 - p][i - 1] == 0) //pivot is zero
        {
            bool P = true;
            for (j = i - p; j < m; ++j) //row exchange is allowed
                if (A[j][i - 1] != 0) {
                    for (l = 0; l < n; ++l)  A[i - 1 - p][l] ^= A[j][l] ^= A[i - 1 - p][l] ^= A[j][l]; //row exchange
                    i--; //row exchanged. do it again
                    P = false;
                    break;
                }
            if (!P)
                continue;
            p++;
            continue;
        }
        piv.push_back(i - 1); //pivot location tracker
        rank++;
        long long temp = A[i - 1 - p][i - 1];
        A[i - 1 - p][i - 1] = 1;
        for (j = i; j < n; ++j)  A[i - 1 - p][j] = (A[i - 1 - p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD - A[j][i - 1];
            if (mul == 0)
                continue;
            for (k = i - 1; k < n; ++k)
                A[j][k] = (A[j][k] + A[i - 1 - p][k] * mul) % MOD;
        }
    }
    if (rank == n) {
        vector<vector<long long>> NR;
        return NR;
    }
    for (i = (int)piv.size() - 1; i > 0; --i) //upper elimination
        for (j = i - 1; j >= 0; --j)
        {
            long long mul = MOD - A[j][piv[i]];
            for (k = piv[i]; k < n; ++k)
                A[j][k] = (A[j][k] + A[i][k] * mul) % MOD;
        }
    for (i = m - 1; i >= 0; --i) //zero row pop
    {
        bool P = false;
        for (j = 0; j < n; ++j)
            if (A[i][j] != 0) {
                P = true;
                break;
            }
        if (P)
            break;
        A.pop_back();
    }
    if (A.empty())
        return I_n(n);
    vector<vector<long long>> TR = matrix_transpose(A), F(n - rank, vector<long long>(rank, 0));
    vector<pair<int, int>> exc;
    for (i = 0; i < rank; ++i) {
        if (piv[i] != i) {
            for (j = 0; j < TR[i].size(); ++j)
                TR[i][j] ^= TR[piv[i]][j] ^= TR[i][j] ^= TR[piv[i]][j];
            exc.push_back({ i,piv[i] });
            piv[i] = i;
        }
    }
    for (p = i; i < n; ++i)
        for (j = 0; j < rank; ++j)
            F[i - p][j] = TR[i][j] ? MOD - TR[i][j] : 0; //remaining col of A is not from I
    vector<vector<long long>> N = matrix_transpose(F);
    for (i = 0; i < F.size(); ++i) {
        vector<long long> te(F.size(), 0);
        te[i] = 1;
        N.push_back(te); // I padding
    }
    for (i = (int)exc.size() - 1; i >= 0; --i)
        for (j = 0; j < N.front().size(); ++j)
            N[exc[i].first][j] ^= N[exc[i].second][j] ^= N[exc[i].first][j] ^= N[exc[i].second][j];
    if (!Orth)
        return N;
    vector<vector<long long>> X = matrix_transpose(N);  //G-S process
    vector<vector<long long>> V = X;
    vector<long long> DP(2, 0);
    DP[1] = (V[0] * V[0]) % MOD;
    for (i = 1; i < N[0].size(); ++i)
    {
        for (j = 1; j <= i; ++j)
        {
            long long c = (V[j - 1] * X[i] % MOD) * inverse(DP[j]) % MOD;
            for (l = 0; l < V[i].size(); ++l)
                V[i][l] = (V[i][l] + MOD - (c * V[j - 1][l] % MOD)) % MOD;
        }
        DP.push_back((V[i] * V[i]) % MOD);
        if (DP.back() == 0) {
            printf("NullSpace's G-S Process Error : Matrix has dependent column\n\n");
            exit(1);
        }
    }
    return matrix_transpose(V);
}
inline vector<long long> Ax_b(vector<vector<long long>>& A, vector<long long> b) {
    if (A.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }
    int m = (int)A.size(), n = (int)A[0].size();
    n++;
    vector<vector<long long>> R(m, vector<long long>(n, 0));
    vector<int> piv;
    int i, j, k, p = 0, elimi = 0;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n - 1; ++j)
            R[i][j] = A[i][j];
        R[i][j] = b[i];
    }
    for (i = 1; i < n && i - 1 - p < m; ++i)
    {
        if (R[i - 1 - p][i - 1] == 0)
        {
            bool P = true;
            for (j = i - p; j < m; ++j)
                if (R[j][i - 1] != 0) {
                    vector<long long> temp = R[i - 1 - p];
                    R[i - 1 - p] = R[j];
                    R[j] = temp;
                    i--;
                    P = false;
                    break;
                }
            if (!P)  continue;
            p++;
            continue;
        }
        piv.push_back(i - 1);
        elimi++;
        long long temp = R[i - 1 - p][i - 1];
        R[i - 1 - p][i - 1] = 1;
        for (j = i; j < n; ++j)    R[i - 1 - p][j] = (R[i - 1 - p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD - R[j][i - 1];
            if (mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                R[j][k] = (R[j][k] + R[i - 1 - p][k] * mul) % MOD;
        }
    }
    for (i = (int)piv.size() - 1; i > 0; --i) //upper elimination
        for (j = i - 1; j >= 0; --j)
        {
            long long mul = MOD - R[j][piv[i]];
            for (k = piv[i]; k < n; ++k)
                R[j][k] = (R[j][k] + R[i][k] * mul) % MOD;
        }
    for (i = m - 1; i >= 0; --i) //zero row solvablity
    {
        bool P = false;
        for (j = 0; j < n - 1; ++j)
            if (R[i][j] != 0) {
                P = true;
                break;
            }
        if (P)   break;
        if (R[i][n - 1] != 0) {
            printf("Ax=b calculation Error : This System is Not Solvable\n\n");
            //exit(1);
            return { -1 };
        }
        R.pop_back();
        b.pop_back();
    }
    if (R.size() == R[0].size() - 1) {
        vector<long long> r(R.size());
        for (i = 0; i < r.size(); ++i)   r[i] = R[i][n - 1];
        return r;
    }
    vector<long long> r(n - 1, 0);
    for (i = 0, p = 0; i < n - 1; ++i) {
        if (piv[p] == i) {
            r[i] = R[p][n - 1];
            p++;
        }
    }
    return r;
}
inline bool is_in(vector<vector<long long>>& A, vector<long long> b) {
    if (A.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }
    int m = (int)A.size(), n = (int)A[0].size();
    n++;
    vector<vector<long long>> R(m, vector<long long>(n, 0));
    vector<int> piv;
    int i, j, k, p = 0, elimi = 0;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n - 1; ++j)
            R[i][j] = A[i][j];
        R[i][j] = b[i];
    }
    for (i = 1; i < n && i - 1 - p < m; ++i)
    {
        if (R[i - 1 - p][i - 1] == 0)
        {
            bool P = true;
            for (j = i - p; j < m; ++j)
                if (R[j][i - 1] != 0) {
                    vector<long long> temp = R[i - 1 - p];
                    R[i - 1 - p] = R[j];
                    R[j] = temp;
                    i--;
                    P = false;
                    break;
                }
            if (!P)  continue;
            p++;
            continue;
        }
        piv.push_back(i - 1);
        elimi++;
        long long temp = R[i - 1 - p][i - 1];
        R[i - 1 - p][i - 1] = 1;
        for (j = i; j < n; ++j)    R[i - 1 - p][j] = (R[i - 1 - p][j] * inverse(temp)) % MOD;
        for (j = i - p; j < m; ++j) {
            long long mul = MOD - R[j][i - 1];
            if (mul == 0)  continue;
            for (k = i - 1; k < n; ++k)
                R[j][k] = (R[j][k] + R[i - 1 - p][k] * mul) % MOD;
        }
    }
    for (i = (int)piv.size() - 1; i > 0; --i) //upper elimination
        for (j = i - 1; j >= 0; --j)
        {
            long long mul = MOD - R[j][piv[i]];
            for (k = piv[i]; k < n; ++k)
                R[j][k] = (R[j][k] + R[i][k] * mul) % MOD;
        }
    for (i = m - 1; i >= 0; --i) //zero row solvablity
    {
        bool P = false;
        for (j = 0; j < n - 1; ++j)
            if (R[i][j] != 0) {
                P = true;
                break;
            }
        if (P)   break;
        if (R[i][n - 1] != 0)
            return false;
        R.pop_back();
        b.pop_back();
    }
    return true;
}

inline void Shamir_Shares_Generate(long long K, long long real, long long fake, long long s, long long t, vector<pair<long long, long long>>& Points, vector<vector<long long>>& M) {
    /*
    1. initiate Shamir-SS for a polynomial degree of n with random secret key K.
    2. create real share and fake shares
    3. create DHF matrix M
    4. Mix positions
    */

    //1
    vector<long long> Shamir_Poly(t+1);
    Shamir_Poly[0] = K;
    for(int i=1; i<Shamir_Poly.size(); ++i)
        Shamir_Poly[i] = rand() % MOD;

    //2
    Points.clear();
    set<long long> x_cor;
    while(x_cor.size() < real+fake)
        x_cor.insert(rand() % MOD);
    for(auto a : x_cor) {
        long long r = K;
        for(int i=1; i<t+1; ++i)
            r = (r + Shamir_Poly[i] * power(a,i)) % MOD;
        Points.push_back({a,r});
    }
    for(int i=0; i<fake; ++i) {
        long long rk = rand() % MOD;
        while(Points[i].second == rk)
            rk = rand() % MOD;
        Points[i].second = rand() % MOD;
    }

    //3
    vector<vector<long long>> KM(s, vector<long long>(t));
    M.clear();
    for(int i=0; i<s; ++i)
        for(int j=0; j<t; ++j)
            KM[i][j] = rand() % MOD;
    for(int i=0; i<fake; ++i) { //fake DHF value
        vector<long long> r(s+1);
        for(int j=0; j<=s; ++j)
            r[j] = rand() % MOD;
        M.push_back(r);
    }
    for(int i=fake; i<real+fake; ++i) { //real DHF value
        vector<long long> r(s+1);
        r[0] = Points[i].first;
        for(long long j=0, temp=KM[j][0]; j<s; ++j) {
            for(int l=1; l<t; ++l)
                temp = (temp + KM[j][l] * power(r[0], l)) % MOD;
            r[j+1] = temp;
        }
        M.push_back(r);
    }

    //4
    for(int i=0; i<real + fake; ++i) {
        int k = rand() % (fake+real);
        auto T1 = Points[i];
        Points[i] = Points[k];
        Points[k] = T1;

        auto T2 = M[i];
        M[i] = M[k];
        M[k] = T2;
    }
}

inline long long Detection_Algorithm(long long s, long long t, vector<pair<long long, long long>>& Points, vector<vector<long long>>& v) {
    vector<vector<long long>> M(s+t, vector<long long>(v.size()));
    long long m = v.size();
    for(int i=0; i<t; ++i)
        for(int j=0; j<m; ++j)
            M[i][j] = power(v[j][0],i);
    for(int i=t; i<t+s; ++i)
        for(int j=0; j<m; ++j)
            M[i][j] = v[j][i-t+1];
    vector<vector<long long>> NL = Null_Space(M, false);
    if(NL.empty())
        return -1;
    vector<long long> b;
    vector<vector<long long>> A;
    for(int i=0; i<m && b.size() <= t+1; ++i)
        if(NL[i][0] != 0) {
            b.push_back(Points[i].second);
            vector<long long> tv;
            for(int j=0; j<t+1; ++j)
                tv.push_back(power(Points[i].first, j));
            A.push_back(tv);
        }
    vector<long long> x = Ax_b(A,b);
    return x[0];
}

inline void testing() {
    double avt = 0;
    for (int trial = 1; trial < 1000000000; ++trial) {
        vector<pair<long long, long long>> v;
        vector<vector<long long>> M;
        long long K = rand() % MOD, s = 20, t = 30, real = 31, fake = 5; //only when real > t, it is possible to reconstruct K 
        bool possible = real > t;

        auto start = chrono::high_resolution_clock::now();
        Shamir_Shares_Generate(K, real, fake, s, t, v, M);
        long long reconstructed_K = Detection_Algorithm(s, t, v, M);
        auto end = chrono::high_resolution_clock::now();

        if ((reconstructed_K == -1 && possible) || (possible && K != reconstructed_K)) {
            printf("NOT GOOD...\n\n");
            printf("original K      = %lld\nrecunstructed K = %lld\n", K, reconstructed_K);

            exit(1);
        }
        printf("original K      = %lld\nrecunstructed K = %lld\n", K, reconstructed_K);
        chrono::duration<double> e1 = end - start;
        double d1 = (double)(e1.count());
        avt += d1;
        printf("-- %d\t\t%lf sec.\t\t(avg %lf sec)\n", trial, d1, avt / trial);
    }
}

int main()
{
    //MOD = 1000000007;         //2*500000003         worst distributed
    //MOD = 100000007;          //2*491*101833
    //MOD = 131071;             //2*3*5*17*257
    MOD = 524287;             //2*3*3*3*7*19*73     well distributed
    //MOD = 65537;              //2^16
    //MOD = 653659;               //2*3*108943
    //MOD = 101;                //2*2*5*5
    Initiation();

    testing();

    vector<pair<long long, long long>> v;
    vector<vector<long long>> M;
    long long K = rand() % MOD, s = 5, t = 5, real = 6, fake = 5; //only when real > t, it is possible to reconstruct K 
    bool possible = real > t;

    auto start = chrono::high_resolution_clock::now();
    Shamir_Shares_Generate(K, real, fake, s, t, v, M);
    long long reconstructed_K = Detection_Algorithm(s, t, v, M);
    auto end = chrono::high_resolution_clock::now();

    if ((reconstructed_K == -1 && possible) || (possible && K != reconstructed_K)) {
        printf("NOT GOOD...\n\n");
        printf("original K      = %lld\nrecunstructed K = %lld\n", K, reconstructed_K);

        exit(1);
    }
    printf("original K      = %lld\nrecunstructed K = %lld\n", K, reconstructed_K);
    chrono::duration<double> e1 = end - start;
    double d1 = (double)(e1.count());
    printf("-- \t\t%lf sec\n", d1);

    return 0;
}
