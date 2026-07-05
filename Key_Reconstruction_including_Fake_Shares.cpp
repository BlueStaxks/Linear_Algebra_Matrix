#include <fstream>
#include <time.h>
#include <chrono>
#include <random>

#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

#include <set>
#include <math.h>
#include <vector>
#include <numeric>
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

mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
inline long long get_rand(long long mod) {
    uniform_int_distribution<long long> dist(0, mod - 1);
    return dist(rng);
}


std::unique_ptr<tbb::spin_mutex[]> inverse_locks(new tbb::spin_mutex[1024]);
inline long long inverse(long long a) {
    if (!a) {
        printf("Integer Inverse Error : 0 has no inverse.\n\n");
        exit(1);
    }
    int lock_index = a % 1024;

    {
        tbb::spin_mutex::scoped_lock lock(inverse_locks[lock_index]);
        if (int_inverse[a]) {
            return int_inverse[a];
        }
    }

    long long q, r1 = MOD, r2 = a, r = 1, t1 = 0, t2 = 1, t;
    while (r) {
        q = r1 / r2;
        r = r1 % r2;
        t = (t1 - q * t2) % MOD;
        if (t < 0) t += MOD;
        r1 = r2;  r2 = r;   t1 = t2;  t2 = t;
    }

    {
        tbb::spin_mutex::scoped_lock lock(inverse_locks[lock_index]);
        int_inverse[a] = (int)t1;
    }

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
inline vector<long long> decompose(long long a) {
    vector<long long> r;
    while (!(a & 1)) {
        r.push_back(2);
        a >>= 1;
    }
    for (long long i = 3; i * i <= a; i += 2) {
        while (a % i == 0) {
            r.push_back(i);
            a /= i;
        }
    }
    if (a > 1)  r.push_back(a);
    return r;
}
inline vector<long long> divisor(long long a) {
    vector<long long> r;
    for (long long i = 1; i * i <= a; ++i)
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
    string name = to_string(MOD) + ".bin";
    ifstream file1(name);

    if (!file1) { // No file
        int_inverse.resize(MOD, 0);
        MOD_decompose = decompose(MOD - 1);
        MOD_divisors = divisor(MOD - 1);

        primitive = 0;
        for (int i = 2; i < MOD; ++i) {
            bool is_primitive = true;
            for (int j = (int)MOD_divisors.size() - 2; j >= 1; --j) {
                if (power(i, MOD_divisors[j]) == 1) {
                    is_primitive = false;
                    break;
                }
            }
            if (is_primitive) {
                primitive = i; // Found the smallest primitive root
                break;
            }
        }

        seeds.resize(MOD);
        ones_roots.resize(MOD);
        for (long long i = 1, t = primitive; i < MOD; ++i, t = (t * primitive) % MOD) {
            seeds[t] = (int)i;  // primitive ^ seeds[i] = i
        }

        // Create an array of spin mutexes. One lock for each divisor index to prevent massive contention.
        std::unique_ptr<tbb::spin_mutex[]> locks(new tbb::spin_mutex[MOD_divisors.size()]);
        tbb::parallel_for(1LL, MOD, [&](long long i) {
            for (int j = 0; j < (int)MOD_divisors.size() - 1; ++j) {
                if (power(i, MOD_divisors[j]) == 1) {
                    // Lock ONLY the specific vector we are writing to
                    tbb::spin_mutex::scoped_lock lock(locks[j]);
                    ones_roots[MOD_divisors[j]].push_back((int)i);
                }
            }
            });

        writeData(name);
    }
    else {
        readData(name);
    }
    printf("Init done.\n");
}


inline void matrix_print(const vector<vector<long long>> a) {
    if (a.empty())   return;   if (a[0].empty())    return;
    for (int i = 0; i < a.size(); ++i)
    {
        for (int j = 0; j < a[0].size(); ++j)
            printf("%lld\t", a[i][j]);
        printf("\n");
    }
    printf("\n\n");
}
inline void vector_print(const vector<long long> a) {
    for (int i = 0; i < a.size(); ++i)
        printf("%lld\t", a[i]);
    printf("\n\n\n");
}


inline vector<vector<long long>> operator * (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty()) {
        return {};
    }
    if (a[0].size() != b.size()) {
        printf("Matrix Multiplication Error : Matrix size does not match\n");
        exit(1);
    }

    int rowsA = a.size();
    int colsA = a[0].size();
    int colsB = b[0].size();

    // Transpose B to prevent Cache Misses
    vector<vector<long long>> b_T(colsB, vector<long long>(colsA));
    for (int i = 0; i < colsB; ++i) {
        for (int j = 0; j < colsA; ++j) {
            b_T[i][j] = b[j][i];
        }
    }

    vector<vector<long long>> R(rowsA, vector<long long>(colsB, 0));
    tbb::parallel_for(0, rowsA, [&](int i) {
        for (int j = 0; j < colsB; ++j) {
            unsigned long long sum = 0;
            for (int k = 0; k < colsA; ++k) {
                sum += (unsigned long long)a[i][k] * b_T[j][k];
                if (k & 1)  sum %= MOD;
            }
            R[i][j] = sum % MOD;
        }
        });

    return R;
}
inline vector<long long> operator * (const vector<vector<long long>>& a, const vector<long long>& b) {
    if (a.empty() || a[0].empty() || b.empty()) {
        return {};
    }
    if (a[0].size() != b.size()) {
        printf("Matrix Vector Multiplication Error : Matrix and Vector's size do not match\n");
        exit(1);
    }

    int rowsA = a.size();
    int colsA = a[0].size();

    vector<long long> R(rowsA, 0);
    tbb::parallel_for(0, rowsA, [&](int i) {
        unsigned long long sum = 0;
        for (int j = 0; j < colsA; ++j) {
            sum += (unsigned long long)a[i][j] * b[j];
            if (j & 1)  sum %= MOD;
        }
        R[i] = sum % MOD;
        });

    return R;
}
inline long long operator * (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Dot Product Error : Vector size does not match\n");
        exit(1);
    }
    unsigned long long sum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += (unsigned long long)a[i] * b[i];
        if (i & 1)  sum %= MOD;
    }
    return sum % MOD;
}
inline vector<long long> operator * (const long long a, const vector<long long>& b) {
    if (a == 0) return vector<long long>(b.size(), 0);
    if (a == 1) return b;
    vector<long long> R(b.size());
    for (auto i = 0; i < b.size(); ++i)
        R[i] = (a * b[i]) % MOD;
    return R;
}
inline vector<vector<long long>> operator * (const long long a, const vector<vector<long long>>& b) {
    int rows = b.size();
    int cols = b[0].size();
    if (a == 0) return vector<vector<long long>>(rows, vector<long long>(cols, 0));
    if (a == 1) return b;
    vector<vector<long long>> R(rows, vector<long long>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            R[i][j] = (a * b[i][j]) % MOD;
    return R;
}
inline vector<vector<long long>> operator + (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty()) return {};
    if (a[0].size() != b[0].size() || a.size() != b.size()) {
        printf("Matrix Addition Error : Matrix size does not match\n");
        exit(1);
    }
    int rows = a.size();
    int cols = a[0].size();
    vector<vector<long long>> R(rows, vector<long long>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            long long val = a[i][j] + b[i][j];
            if (val >= MOD) val -= MOD;
            R[i][j] = val;
        }
    }
    return R;
}
inline vector<vector<long long>> operator - (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty()) return {};
    if (a[0].size() != b[0].size() || a.size() != b.size()) {
        printf("Matrix Subtraction Error : Matrix size does not match\n");
        exit(1);
    }
    int rows = a.size();
    int cols = a[0].size();
    vector<vector<long long>> R(rows, vector<long long>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            long long val = a[i][j] - b[i][j];
            if (val < 0) val += MOD;
            R[i][j] = val;
        }
    }
    return R;
}
inline vector<vector<long long>> operator | (const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || a[0].empty()) return b;
    if (b.empty() || b[0].empty()) return a;

    int rowsA = a.size();
    int colsA = a[0].size();
    int rowsB = b.size();
    int colsB = b[0].size();

    vector<vector<long long>> R(rowsA + rowsB, vector<long long>(colsA + colsB, 0));
    for (int i = 0; i < rowsA; ++i)
        copy(a[i].begin(), a[i].end(), R[i].begin());
    for (int i = 0; i < rowsB; ++i)
        copy(b[i].begin(), b[i].end(), R[rowsA + i].begin() + colsA);
    return R;
}
inline vector<long long> operator + (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Addition Error : Vector size does not match\n");
        exit(1);
    }
    vector<long long> R(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        long long val = a[i] + b[i];
        if (val >= MOD) val -= MOD;
        R[i] = val;
    }
    return R;
}
inline vector<long long> operator - (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        printf("Vector Subtraction Error : Vector size does not match\n");
        exit(1);
    }
    vector<long long> R(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        long long val = a[i] - b[i];
        if (val < 0) val += MOD;
        R[i] = val;
    }
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
    for (size_t i = 0; i < n; ++i)  I[i][i] = 1;
    return I;
}
inline vector<vector<long long>> I_n(int n, long long a) {
    vector<vector<long long>> I(n, vector<long long>(n, 0));
    for (size_t i = 0; i < n; ++i)  I[i][i] = a;
    return I;
}
inline vector<vector<long long>> matrix_transpose(const vector<vector<long long>>& a) {
    if (a.empty() || a[0].empty()) return {};
    vector<vector<long long>> R(a.front().size(), vector<long long>(a.size()));
    size_t i, j;
    for (i = 0; i < a.size(); ++i)
        for (j = 0; j < a.front().size(); ++j)
            R[j][i] = a[i][j];
    return R;
}
inline vector<vector<long long>> matrix_transpose_tiled(const vector<vector<long long>>& a) {
    if (a.empty() || a[0].empty()) return {};
    size_t rows = a.size();
    size_t cols = a[0].size();
    vector<vector<long long>> R(cols, vector<long long>(rows));
    const size_t BLOCK_SIZE = 32;
    for (size_t i = 0; i < rows; i += BLOCK_SIZE) {
        for (size_t j = 0; j < cols; j += BLOCK_SIZE) {
            size_t max_i = min(i + BLOCK_SIZE, rows);
            size_t max_j = min(j + BLOCK_SIZE, cols);

            for (size_t ii = i; ii < max_i; ++ii)
                for (size_t jj = j; jj < max_j; ++jj)
                    R[jj][ii] = a[ii][jj];
        }
    }
    return R;
}
inline vector<vector<long long>> matrix_power(vector<vector<long long>> a, unsigned long long n) {
    vector<vector<long long>> res = I_n(a.size());
    while (n) {
        if (n & 1)  res = res * a;
        n >>= 1;
        if (!n) break;
        a = a * a;
    }
    return res;
}
inline void matrix_chop(vector<vector<vector<long long>>>& M, const vector<vector<long long>>& F, const vector<int>& list) {
    M.reserve(M.size() + list.size());
    size_t p = 0;
    for (size_t i = 0; i < list.size(); ++i) {
        size_t block_size = list[i];
        if (p + block_size > F.size() || p + block_size > F[0].size()) {
            printf("Matrix Chop Error : Block bounds exceed matrix dimensions\n");
            exit(1);
        }

        vector<vector<long long>> block(block_size, vector<long long>(block_size));
        for (size_t j = 0; j < block_size; ++j) {
            // Source Start: F[row].begin() + column_offset
            // Source End: F[row].begin() + column_offset + block_size
            // Destination: block[j].begin()
            std::copy(
                F[j + p].begin() + p,
                F[j + p].begin() + p + block_size,
                block[j].begin()
            );
        }
        M.push_back(std::move(block));
        p += block_size;
    }
}
inline vector<vector<long long>> matrix_partial_multiply(const vector<vector<long long>>& A, const vector<vector<long long>>& B, const vector<int>& list) {
    size_t n = A.size();
    vector<vector<long long>> R(n, vector<long long>(n, 0));
    size_t p = 0;

    for (size_t i = 0; i < list.size(); ++i) {
        size_t block_size = list[i];
        if (p + block_size > n || p + block_size > B[0].size()) {
            printf("Matrix Partial Multiply Error : Block bounds exceed matrix dimensions\n");
            exit(1);
        }

        vector<vector<long long>> B_block_T(block_size, vector<long long>(block_size));
        for (size_t r = 0; r < block_size; ++r) {
            for (size_t c = 0; c < block_size; ++c) {
                B_block_T[c][r] = B[p + r][p + c];
            }
        }

        for (size_t j = 0; j < n; ++j) {
            for (size_t l = 0; l < block_size; ++l) {
                unsigned long long sum = 0;
                for (size_t k = 0; k < block_size; ++k) {
                    sum += (unsigned long long)A[j][p + k] * B_block_T[l][k];
                    if (k & 1)  sum %= MOD;
                }
                R[j][p + l] = (R[j][p + l] + sum) % MOD;
            }
        }
        p += block_size;
    }
    return R;
}
inline size_t matrix_rank(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D[0].empty()) return 0;

    size_t m = A_2D.size();
    size_t n = A_2D[0].size();

    vector<long long> A(m * n);
    for (size_t i = 0; i < m; ++i)
        copy(A_2D[i].begin(), A_2D[i].end(), A.begin() + i * n);

    size_t rank = 0;
    size_t row = 0; // Tracks the current pivot row
    for (size_t col = 0; col < n && row < m; ++col) {
        size_t pivot_row = row;
        while (pivot_row < m && A[pivot_row * n + col] == 0) pivot_row++;
        if (pivot_row == m) continue; // All zeros in this column below 'row'
        if (pivot_row != row)
            swap_ranges(A.begin() + row * n, A.begin() + (row + 1) * n, A.begin() + pivot_row * n);

        long long inv = inverse(A[row * n + col]);
        for (size_t k = col; k < n; ++k)
            A[row * n + k] = (A[row * n + k] * inv) % MOD;

        tbb::parallel_for(tbb::blocked_range<size_t>(row + 1, m),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t j = r.begin(); j != r.end(); ++j) {
                    long long factor = A[j * n + col];
                    if (factor == 0) continue; // Optimization: Skip if already 0
                    for (size_t k = col; k < n; ++k) {
                        long long sub = (A[row * n + k] * factor) % MOD;
                        A[j * n + k] -= sub;
                        if (A[j * n + k] < 0) A[j * n + k] += MOD;
                    }
                }
            });
        row++;
        rank++;
    }
    return rank;
}
inline vector<vector<long long>> matrix_inverse(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D.size() != A_2D.front().size()) {
        printf("Matrix Inversion Error : Matrix is not square\n\n");
        exit(1);
    }
    size_t n = A_2D.size();
    vector<long long> A(n * n);
    vector<long long> I(n * n, 0);
    for (size_t i = 0; i < n; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), A.begin() + i * n);
        I[i * n + i] = 1; // Initialize Identity matrix concurrently
    }

    for (size_t p = 0; p < n - 1; ++p) {
        if (A[p * n + p] == 0) {
            bool found = false;
            for (size_t j = p + 1; j < n; ++j) {
                if (A[j * n + p] != 0) {
                    swap_ranges(A.begin() + j * n, A.begin() + (j + 1) * n, A.begin() + p * n);
                    swap_ranges(I.begin() + j * n, I.begin() + (j + 1) * n, I.begin() + p * n);
                    found = true;
                    break;
                }
            }
            if (!found) {
                printf("Matrix Inversion Error : Matrix is singular\n\n");
                exit(1);
            }
        }

        long long inv_pivot = inverse(A[p * n + p]);
        tbb::parallel_for(tbb::blocked_range<size_t>(p + 1, n),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t j = r.begin(); j != r.end(); ++j) {
                    long long mul = (MOD - A[j * n + p]) * inv_pivot % MOD;
                    for (size_t k = p; k < n; ++k)    // Skip zeros before 'p'
                        A[j * n + k] = (A[j * n + k] + A[p * n + k] * mul) % MOD;
                    for (size_t k = 0; k < n; ++k)
                        I[j * n + k] = (I[j * n + k] + I[p * n + k] * mul) % MOD;
                }
            });
    }

    if (n > 0 && A[(n - 1) * n + (n - 1)] == 0) {
        printf("Matrix Inversion Error : Matrix is singular\n\n");
        exit(1);
    }
    for (int p = n - 1; p > 0; --p) {
        long long inv_pivot = inverse(A[p * n + p]);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, p),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t j = r.begin(); j != r.end(); ++j) {
                    long long mul = (MOD - A[j * n + p]) * inv_pivot % MOD;
                    A[j * n + p] = (A[j * n + p] + A[p * n + p] * mul) % MOD;   // Skip zeros entirely, only process column 'p'
                    for (size_t k = 0; k < n; ++k)
                        I[j * n + k] = (I[j * n + k] + I[p * n + k] * mul) % MOD;
                }
            });
    }

    vector<vector<long long>> I_out(n, vector<long long>(n));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                long long t = inverse(A[i * n + i]);
                for (size_t j = 0; j < n; ++j)
                    I_out[i][j] = (I[i * n + j] * t) % MOD;
            }
        });
    return I_out;
}
inline long long matrix_determinant(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D.size() != A_2D.front().size()) {
        printf("Matrix determinant Error : Matrix is not square\n\n");
        exit(1);
    }
    size_t n = A_2D.size();
    vector<long long> A(n * n);
    for (size_t i = 0; i < n; ++i)
        copy(A_2D[i].begin(), A_2D[i].end(), A.begin() + i * n);
    long long det = 1;

    for (size_t p = 0; p < n; ++p) {
        if (A[p * n + p] == 0) {
            bool found = false;
            for (size_t j = p + 1; j < n; ++j) {
                if (A[j * n + p] != 0) {
                    swap_ranges(A.begin() + j * n, A.begin() + (j + 1) * n, A.begin() + p * n);
                    det = MOD - det;
                    found = true;
                    break;
                }
            }
            if (!found) return 0;
        }

        det = (det * A[p * n + p]) % MOD;
        if (p == n - 1) break;

        long long inv_pivot = inverse(A[p * n + p]);
        tbb::parallel_for(tbb::blocked_range<size_t>(p + 1, n),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t j = r.begin(); j != r.end(); ++j) {
                    long long factor = A[j * n + p];
                    if (factor == 0) continue; // Skip already zeroed rows
                    long long mul = (MOD - factor) * inv_pivot % MOD;
                    for (size_t k = p + 1; k < n; ++k)
                        A[j * n + k] = (A[j * n + k] + A[p * n + k] * mul) % MOD;
                }
            });
    }
    return det;
}
inline vector<vector<long long>> Null_Space(const vector<vector<long long>>& A_2D, bool Orth) {
    if (A_2D.empty() || A_2D[0].empty()) return {};

    int m = A_2D.size();
    int n = A_2D[0].size();
    vector<long long> A(m * n);
    for (int i = 0; i < m; ++i)
        copy(A_2D[i].begin(), A_2D[i].end(), A.begin() + i * n);
    vector<int> piv;
    int row = 0;
    for (int col = 0; col < n && row < m; ++col) {
        int pivot_row = row;
        while (pivot_row < m && A[pivot_row * n + col] == 0) pivot_row++;
        if (pivot_row == m) continue; // Free variable column
        if (pivot_row != row)
            swap_ranges(A.begin() + row * n, A.begin() + (row + 1) * n, A.begin() + pivot_row * n);

        piv.push_back(col);
        long long inv = inverse(A[row * n + col]);
        for (int k = col; k < n; ++k)
            A[row * n + k] = (A[row * n + k] * inv) % MOD;

        tbb::parallel_for(tbb::blocked_range<int>(row + 1, m),
            [&](const tbb::blocked_range<int>& r) {
                for (int j = r.begin(); j != r.end(); ++j) {
                    long long factor = A[j * n + col];
                    if (factor == 0) continue;
                    for (int k = col; k < n; ++k) {
                        long long sub = (A[row * n + k] * factor) % MOD;
                        A[j * n + k] -= sub;
                        if (A[j * n + k] < 0)   A[j * n + k] += MOD;
                    }
                }
            });
        row++;
    }

    int rank = row;
    if (rank == n) return {}; // Trivial null space (empty)

    for (int i = rank - 1; i >= 0; --i) {
        int p_col = piv[i];
        tbb::parallel_for(tbb::blocked_range<int>(0, i),
            [&](const tbb::blocked_range<int>& r) {
                for (int j = r.begin(); j != r.end(); ++j) {
                    long long factor = A[j * n + p_col];
                    if (factor == 0) continue;
                    for (int k = p_col; k < n; ++k) {
                        long long sub = (A[i * n + k] * factor) % MOD;
                        A[j * n + k] -= sub;
                        if (A[j * n + k] < 0)   A[j * n + k] += MOD;
                    }
                }
            });
    }

    int null_dim = n - rank;
    vector<vector<long long>> NS(null_dim, vector<long long>(n, 0));

    vector<bool> is_pivot(n, false);
    for (int p : piv) is_pivot[p] = true;

    int free_idx = 0;
    for (int col = 0; col < n; ++col) {
        if (!is_pivot[col]) {
            NS[free_idx][col] = 1; // The identity portion
            for (int i = 0; i < rank; ++i)
                NS[free_idx][piv[i]] = (MOD - A[i * n + col]) % MOD;
            free_idx++;
        }
    }

    if (Orth) {
        vector<long long> DP(null_dim, 0);
        for (int i = 0; i < null_dim; ++i) {
            for (int j = 0; j < i; ++j) {
                long long dot = 0;
                for (int k = 0; k < n; ++k)
                    dot = (dot + NS[i][k] * NS[j][k]) % MOD;

                long long c = (dot * inverse(DP[j])) % MOD;

                for (int k = 0; k < n; ++k) {
                    long long sub = (c * NS[j][k]) % MOD;
                    NS[i][k] -= sub;
					if (NS[i][k] < 0)   NS[i][k] += MOD;
                }
            }
            long long norm_sq = 0;
            for (int k = 0; k < n; ++k)
                norm_sq = (norm_sq + NS[i][k] * NS[i][k]) % MOD;
            DP[i] = norm_sq;
            if (DP[i] == 0) {
                printf("NullSpace's G-S Process Error : Isotropic vector encountered (v*v = 0 mod P)\n\n");
                exit(1);
            }
        }
    }

    vector<vector<long long>> NS_col(n, vector<long long>(null_dim));
    tbb::parallel_for(tbb::blocked_range<int>(0, null_dim),
        [&](const tbb::blocked_range<int>& r) {
            for (int i = r.begin(); i != r.end(); ++i)
                for (int j = 0; j < n; ++j)
                    NS_col[j][i] = NS[i][j];
        });

    return NS_col;
}
inline vector<long long> Ax_b(const vector<vector<long long>>& A_2D, const vector<long long>& b) {
    if (A_2D.empty() || A_2D.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }

    int m = A_2D.size();
    int n = A_2D[0].size();
    int cols = n + 1; // Width of the augmented matrix [A | b]

    vector<long long> R(m * cols, 0);
    for (int i = 0; i < m; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), R.begin() + i * cols);
        R[i * cols + n] = b[i];
    }
    vector<int> piv;
    int row = 0;
    for (int col = 0; col < n && row < m; ++col) {
        int pivot_row = row;
        while (pivot_row < m && R[pivot_row * cols + col] == 0) pivot_row++;
        if (pivot_row == m) continue; // Free variable
        if (pivot_row != row)
            swap_ranges(R.begin() + row * cols, R.begin() + (row + 1) * cols, R.begin() + pivot_row * cols);

        piv.push_back(col);
        long long inv = inverse(R[row * cols + col]);
        for (int k = col; k < cols; ++k)
            R[row * cols + k] = (R[row * cols + k] * inv) % MOD;

        tbb::parallel_for(tbb::blocked_range<int>(row + 1, m),
            [&](const tbb::blocked_range<int>& r) {
                for (int j = r.begin(); j != r.end(); ++j) {
                    long long factor = R[j * cols + col];
                    if (factor == 0) continue;
                    for (int k = col; k < cols; ++k) {
                        long long sub = (R[row * cols + k] * factor) % MOD;
                        R[j * cols + k] -= sub;
						if (R[j * cols + k] < 0) R[j * cols + k] += MOD;
                    }
                }
            });
        row++;
    }

    int rank = row;
    for (int i = rank; i < m; ++i) {
        if (R[i * cols + n] != 0) {
            printf("Ax=b calculation Error : This System is Not Solvable\n\n");
            return vector<long long>(n, 0);
        }
    }

    for (int i = rank - 1; i >= 0; --i) {
        int p_col = piv[i];
        tbb::parallel_for(tbb::blocked_range<int>(0, i),
            [&](const tbb::blocked_range<int>& r) {
                for (int j = r.begin(); j != r.end(); ++j) {
                    long long factor = R[j * cols + p_col];
                    if (factor == 0) continue;

                    // We only need to eliminate the specific column 'p_col' and update 'b'
                    R[j * cols + p_col] -= (R[i * cols + p_col] * factor) % MOD;
                    if (R[j * cols + p_col] < 0) R[j * cols + p_col] += MOD;
                    R[j * cols + n] -= (R[i * cols + n] * factor) % MOD;
                    if (R[j * cols + n] < 0) R[j * cols + n] += MOD;
                }
            });
    }

    vector<long long> x(n, 0);
    for (int i = 0; i < rank; ++i)
        x[piv[i]] = R[i * cols + n];
    return x;
}
inline bool is_in(const vector<vector<long long>>& A_2D, const vector<long long>& b) {
    if (A_2D.empty() || A_2D.size() != b.size()) {
        printf("Ax=b calculation Error : Size is different\n\n");
        exit(1);
    }

    int m = A_2D.size();
    int n = A_2D[0].size();
    int cols = n + 1; // Width of [A | b]

    vector<long long> R(m * cols, 0);
    for (int i = 0; i < m; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), R.begin() + i * cols);
        R[i * cols + n] = b[i]; // Append b
    }
    int row = 0;

    for (int col = 0; col < n && row < m; ++col) {
        int pivot_row = row;
        while (pivot_row < m && R[pivot_row * cols + col] == 0) pivot_row++;
        if (pivot_row == m) continue; // Free variable column
        if (pivot_row != row)
            swap_ranges(R.begin() + row * cols, R.begin() + (row + 1) * cols, R.begin() + pivot_row * cols);

        long long inv = inverse(R[row * cols + col]);
        for (int k = col; k < cols; ++k)
            R[row * cols + k] = (R[row * cols + k] * inv) % MOD;

        tbb::parallel_for(tbb::blocked_range<int>(row + 1, m),
            [&](const tbb::blocked_range<int>& r) {
                for (int j = r.begin(); j != r.end(); ++j) {
                    long long factor = R[j * cols + col];
                    if (factor == 0) continue;
                    for (int k = col; k < cols; ++k) {
                        long long sub = (R[row * cols + k] * factor) % MOD;
                        R[j * cols + k] -= sub;
						if (R[j * cols + k] < 0) R[j * cols + k] += MOD;
                    }
                }
            });
        row++;
    }

    for (int i = row; i < m; ++i)
        if (R[i * cols + n] != 0)
            return false;
    return true;
}
inline void matrix_diagonalize_2x2(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    S = { {0, 0}, {0, 0} };
    D = { {0, 0}, {0, 0} };
    long long a = A[0][0], b = A[0][1];
    long long c = A[1][0], d = A[1][1];
    long long tr = a + d;
    if (tr >= MOD) tr -= MOD;

    long long inv2 = inverse(2);
    long long fr = (tr * inv2) % MOD;
    long long det = (a * d - b * c) % MOD;
    if (det < 0) det += MOD; // Handle negative modulo in C++

    long long inroot = (fr * fr - det) % MOD;
    if (inroot < 0) inroot += MOD;
    if (!inroot)
        D[0][0] = D[1][1] = fr;
    else {
        long long k = seeds[inroot];
        if (k & 1) {        // If the exponent is odd, the square root DOES NOT EXIST in F_p.
            printf("Matrix Diagonalize Error : Discriminant is a non-residue. Eigenvalues exist in F_p^2.\n\n");
            exit(1);
        }
        long long seed2 = power(primitive, k >> 1);
        long long d1 = fr + seed2;
        if (d1 >= MOD) d1 -= MOD;
        long long d2 = fr - seed2;
        if (d2 < 0) d2 += MOD;
        D[0][0] = d1;
        D[1][1] = d2;
    }

    if (D[0][0] == D[1][1]) {
        long long a_minus_l = a - D[0][0];
        if (a_minus_l < 0) a_minus_l += MOD;
        long long d_minus_l = d - D[0][0];
        if (d_minus_l < 0) d_minus_l += MOD;
        if (a_minus_l == 0 && b == 0 && c == 0 && d_minus_l == 0) {
            S[0][0] = 1; S[1][1] = 1; // Identity
        }
        else {
            printf("Matrix Diagonalize Error : 2x2 matrix is defective (not diagonalizable)\n\n");
            exit(1);
        }
    }
    else {
        long long a1 = a - D[0][0];
        if (a1 < 0) a1 += MOD;
        if (a1 != 0 || b != 0) {
            S[0][0] = b;
            S[1][0] = (a1 == 0) ? 0 : MOD - a1; // Fast Negation
        }
        else {
            long long d1 = d - D[0][0];
            if (d1 < 0) d1 += MOD;
            S[0][0] = d1;
            S[1][0] = (c == 0) ? 0 : MOD - c;
        }

        long long a2 = a - D[1][1];
        if (a2 < 0) a2 += MOD;
        if (a2 != 0 || b != 0) {
            S[0][1] = b;
            S[1][1] = (a2 == 0) ? 0 : MOD - a2;
        }
        else {
            long long d2 = d - D[1][1];
            if (d2 < 0) d2 += MOD;
            S[0][1] = d2;
            S[1][1] = (c == 0) ? 0 : MOD - c;
        }
    }
}
inline void matrix_diagonalize_BF(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    int n = A.size();
    if (n == 0 || A.front().size() != n) {
        printf("Matrix diagonalization Error : Matrix is not square\n\n");
        exit(1);
    }

    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    long long trace = 0;
    for (int i = 0; i < n; ++i)
        trace = (trace + A[i][i]) % MOD;

    int vc = 0; // Vector count (number of eigenvectors found so far)
    for (long long lambda = 0; lambda < MOD; ++lambda) {
        vector<vector<long long>> ZN = Null_Space(A, Orth);
        if (!ZN.empty()) {
            int multiplicity = ZN[0].size();
            for (int k = 0; k < multiplicity; ++k) {
                D[vc + k][vc + k] = lambda;
                for (int row = 0; row < n; ++row)
                    S[row][vc + k] = ZN[row][k];
            }
            vc += multiplicity;
            if (vc == n) return; // Fully diagonalized
            if (vc == n - 1) {
                long long eig_sum = 0;
                for (int k = 0; k < n - 1; ++k)
                    eig_sum = (eig_sum + D[k][k]) % MOD;
                long long last_lambda = (trace + MOD - eig_sum) % MOD;
                D[n - 1][n - 1] = last_lambda;
                long long diff = (lambda + MOD - last_lambda) % MOD;
                for (int k = 0; k < n; ++k)
                    A[k][k] = (A[k][k] + diff) % MOD;
                ZN = Null_Space(A, Orth);
                for (int row = 0; row < n; ++row)
                    S[row][n - 1] = ZN[row][0];
                return;
            }
        }
        for (int j = 0; j < n; ++j)
            A[j][j] = (A[j][j] + MOD - 1) % MOD;
    }
    printf("Matrix Diagonalize Error : Matrix is defective (not diagonalizable over F_p)\n\n");
    exit(1);
}
inline void matrix_diagonalize_henry(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    int n = (int)A.size();
    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    vector<vector<long long>> AP_1 = matrix_power(A, MOD - 1);
    for (int i = 0; i < n; ++i)
        AP_1[i][i] = (AP_1[i][i] + MOD - 1) % MOD;
    vector<vector<long long>> ZN1 = Null_Space(AP_1, Orth);
    int eigvec_count = ZN1.empty() ? 0 : (int)ZN1[0].size();
    if (eigvec_count > 0)
        for (int j = 0; j < n; ++j)
            copy(ZN1[j].begin(), ZN1[j].end(), S[j].begin());

    for (int i = 0; i < n; ++i)
        AP_1[i][i] = (AP_1[i][i] + 1) % MOD;
    vector<vector<long long>> ZN2 = Null_Space(AP_1, Orth);
    if (!ZN2.empty()) {
        int zero_count = (int)ZN2[0].size();
        for (int j = 0; j < n; ++j)
            copy(ZN2[j].begin(), ZN2[j].end(), S[j].begin() + eigvec_count);
    }

    vector<vector<long long>> New_A(eigvec_count, vector<long long>(eigvec_count, 0));
    if (eigvec_count > 0 && eigvec_count < n) {
        vector<vector<long long>> S_inv = matrix_inverse(S);
        vector<vector<long long>> AS_left(n, vector<long long>(eigvec_count, 0));

        tbb::parallel_for(tbb::blocked_range<int>(0, n),
            [&](const tbb::blocked_range<int>& r) {
                for (int i = r.begin(); i != r.end(); ++i) {
                    for (int j = 0; j < eigvec_count; ++j) {
                        long long sum = 0;
                        for (int k = 0; k < n; ++k)
                            sum = (sum + A[i][k] * S[k][j]) % MOD;
                        AS_left[i][j] = sum;
                    }
                }
            });
        tbb::parallel_for(tbb::blocked_range<int>(0, eigvec_count),
            [&](const tbb::blocked_range<int>& r) {
                for (int i = r.begin(); i != r.end(); ++i) {
                    for (int j = 0; j < eigvec_count; ++j) {
                        long long sum = 0;
                        for (int k = 0; k < n; ++k)
                            sum = (sum + S_inv[i][k] * AS_left[k][j]) % MOD;
                        New_A[i][j] = sum;
                    }
                }
            });
    }
    else if (eigvec_count == n) {
        New_A = A;
    }

    n = eigvec_count;
    vector<vector<long long>> Ss = I_n(n);
    vector<vector<vector<long long>>> M;
    M.push_back(New_A);
    vector<long long> FE(1, 1);
    long long powC = MOD - 1;

    int mat_i = 0; // Initialize global queue iterator

    for (int pi = 0; pi < MOD_decompose.size(); ++pi) {
        int start_mat_i = mat_i;
        int mati_upperbound = (int)M.size();

        // 1. Pre-calculate 'stp' offsets
        vector<int> stp_offsets(mati_upperbound - start_mat_i, 0);
        int current_stp = 0;
        for (int k = start_mat_i; k < mati_upperbound; ++k) {
            stp_offsets[k - start_mat_i] = current_stp;
            current_stp += M[k].size();
        }

        // Struct to hold results from each parallel block locally
        struct BlockResult {
            vector<vector<vector<long long>>> new_matrices;
            vector<long long> new_FEs;
        };
        vector<BlockResult> block_results(mati_upperbound - start_mat_i);
        powC /= MOD_decompose[pi];
        vector<vector<long long>> ST(n, vector<long long>(n, 0));
        tbb::parallel_for(tbb::blocked_range<int>(start_mat_i, mati_upperbound),
            [&](const tbb::blocked_range<int>& r) {
                for (int m_idx = r.begin(); m_idx != r.end(); ++m_idx) {

                    int local_stp = stp_offsets[m_idx - start_mat_i];
                    auto& local_result = block_results[m_idx - start_mat_i];

                    if (M[m_idx].size() == 1) {
                        local_result.new_matrices.push_back(M[m_idx]);
                        local_result.new_FEs.push_back(M[m_idx][0][0]);
                        ST[local_stp][local_stp] = 1;
                        continue;
                    }
                    if (M[m_idx].size() == 2) {
                        vector<vector<long long>> D2, S2;
                        matrix_diagonalize_2x2(M[m_idx], S2, D2, Orth);
                        local_result.new_matrices.push_back({ {D2[0][0]} });
                        local_result.new_matrices.push_back({ {D2[1][1]} });
                        local_result.new_FEs.push_back(D2[0][0]);
                        local_result.new_FEs.push_back(D2[1][1]);
                        ST[local_stp][local_stp] = S2[0][0];
                        ST[local_stp][local_stp + 1] = S2[0][1];
                        ST[local_stp + 1][local_stp] = S2[1][0];
                        ST[local_stp + 1][local_stp + 1] = S2[1][1];
                        continue;
                    }

                    int m_size = M[m_idx].size();
                    vector<vector<long long>> PM = matrix_power(M[m_idx], powC);
                    long long seed = seeds[FE[m_idx]] * inverse(MOD_decompose[pi]) % MOD;
                    long long seed2 = power(primitive, seed);

                    struct RootResult {
                        long long candidate;
                        vector<vector<long long>> ZN;
                    };
                    tbb::concurrent_vector<RootResult> valid_roots;
                    std::atomic<int> local_eigvec_count{ 0 };

                    int roots_size = ones_roots[MOD_decompose[pi]].size();
                    tbb::parallel_for(tbb::blocked_range<int>(0, roots_size),
                        [&](const tbb::blocked_range<int>& inner_r) {
                            for (int i = inner_r.begin(); i != inner_r.end(); ++i) {
                                if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;

                                long long candidate = seed2 * ones_roots[MOD_decompose[pi]][i] % MOD;
                                vector<vector<long long>> query = PM;
                                for (int j = 0; j < query.size(); ++j) {
                                    query[j][j] -= candidate;
                                    if (query[j][j] < 0) query[j][j] += MOD;
                                }

                                vector<vector<long long>> ZN = Null_Space(query, Orth);

                                if (!ZN.empty()) {
                                    valid_roots.push_back({ candidate, ZN });
                                    local_eigvec_count.fetch_add(ZN[0].size(), std::memory_order_relaxed);
                                }
                            }
                        });

                    vector<vector<long long>> St(m_size, vector<long long>(m_size, 0));
                    int current_col = 0;
                    vector<int> eigspace_dim;

                    for (auto& root : valid_roots) {
                        if (current_col >= m_size) break;

                        local_result.new_FEs.push_back(root.candidate);
                        eigspace_dim.push_back((int)root.ZN[0].size());

                        for (int j = 0; j < root.ZN.size(); ++j)
                            copy(root.ZN[j].begin(), root.ZN[j].end(), St[j].begin() + current_col);
                        current_col += root.ZN[0].size();
                    }

                    vector<vector<long long>> mt = matrix_inverse(St) * M[m_idx] * St;
                    for (int i = 0; i < St.size(); ++i)
                        copy(St[i].begin(), St[i].end(), ST[i + local_stp].begin() + local_stp);

                    matrix_chop(local_result.new_matrices, mt, eigspace_dim);
                }
            }); // END OF tbb::parallel_for

        mat_i = mati_upperbound; // Advance queue iterator to new elements
        for (const auto& res : block_results) {
            for (const auto& mat : res.new_matrices) M.push_back(mat);
            for (auto fe : res.new_FEs) FE.push_back(fe);
        }
        Ss = Ss * ST; // update S sequentially
    }
    for (int Di = 0; mat_i < M.size(); ++mat_i)
        for (int i = 0; i < M[mat_i].size(); ++i, ++Di)
            D[Di][Di] = FE[mat_i];
    vector<vector<long long>> St_final = I_n((int)S.size());
    for (int i = 0; i < Ss.size(); ++i)
        copy(Ss[i].begin(), Ss[i].end(), St_final[i].begin());

    S = S * St_final;
}




inline void func1() {
    int N = 10, i, j, k;
    double avt = 0;
    vector<vector<long long>> I, S1, D1, S2, D2;
    I = I_n(N);
    for (int trial = 1; trial < 1000000000; ++trial) {
        vector<vector<long long>> tm = I;
        for (i = 0; i < N - 1; ++i) {
            for (j = i + 1; j < N; ++j) {
                long long mul = rand() % MOD;
                for (k = 0; k < N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        for (i = N - 1; i > 0; --i) {
            for (j = i - 1; j >= 0; --j) {
                long long mul = rand() % MOD;
                for (k = 0; k < N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        vector<vector<long long>> E(N, vector<long long>(N, 0));
        for (i = 0; i < N; ++i)  E[i][i] = rand() % MOD;
        vector<vector<long long>> DC = tm * E * matrix_inverse(tm);
        clock_t s1, f1, s2, f2;
        s1 = clock();
        matrix_diagonalize_BF(DC, S1, D1, false);
        f1 = clock();
        s2 = clock();
        matrix_diagonalize_henry(DC, S2, D2, false);
        f2 = clock();
        if (DC * S2 != S2 * D2 || matrix_determinant(S2) == 0) {
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
        double d1 = (double)(f1 - s1) / CLOCKS_PER_SEC, d2 = (double)(f2 - s2) / CLOCKS_PER_SEC;
        avt += d1 / d2;
        printf("-- %d\t\t%lf sec vs %lf sec.\t\t%lf time faster!!\t\t(avg %lf time faster)\n", trial, d1, d2, d1 / d2, avt / trial);
        //matrix_print(tm);
    }
}

inline void func2() {
    int N = 600, i, j, k;
    double avt = 0;
    vector<vector<long long>> I, S1, D1, S2, D2;
    I = I_n(N);
    for (int trial = 1; trial < 1000000000; ++trial) {
        vector<vector<long long>> tm = I;
        for (i = 0; i < N - 1; ++i) {
            for (j = i + 1; j < N; ++j) {
                long long mul = rand() % MOD;
                for (k = 0; k < N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        for (i = N - 1; i > 0; --i) {
            for (j = i - 1; j >= 0; --j) {
                long long mul = rand() % MOD;
                for (k = 0; k < N; ++k)
                    tm[j][k] = (tm[j][k] + tm[i][k] * mul) % MOD;
            }
        }
        vector<vector<long long>> E(N, vector<long long>(N, 0));
        for (i = 0; i < N; ++i)  E[i][i] = rand() % MOD;
        vector<vector<long long>> DC = tm * E * matrix_inverse(tm);
        clock_t s2, f2;
        s2 = clock();
        matrix_diagonalize_henry(DC, S2, D2, false);
        f2 = clock();
        if (DC * S2 != S2 * D2 || matrix_determinant(S2) == 0) {
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
        double d2 = (double)(f2 - s2) / CLOCKS_PER_SEC;
        avt += d2;
        printf("-- %d\t\t%lf sec.\t\t(avg %lf sec)\n", trial, d2, avt / trial);
        //matrix_print(tm);
    }
}

inline void Shamir_Shares_Generate(long long K, long long real, long long fake, long long s, long long t, vector<pair<long long, long long>>& Points, vector<vector<long long>>& M) {
    /*
    1. initiate Shamir-SS for a polynomial degree of t with random secret key K.
    2. create real share and fake shares
    3. create DHF matrix M
    4. Mix positions
    */

    // 1
    vector<long long> Shamir_Poly(t + 1);
    Shamir_Poly[0] = K;
    for (int i = 1; i < Shamir_Poly.size(); ++i)
        Shamir_Poly[i] = get_rand(MOD);

    // 2
    Points.clear();
    set<long long> x_cor;
    while (x_cor.size() < real + fake)
        x_cor.insert(get_rand(MOD));

    for (auto a : x_cor) {
        long long r = K;
        for (int i = 1; i < t + 1; ++i)
            r = (r + Shamir_Poly[i] * power(a, i)) % MOD;
        Points.push_back({ a, r });
    }

    for (int i = 0; i < fake; ++i) {
        long long rk = get_rand(MOD);
        while (Points[i].second == rk)
            rk = get_rand(MOD);
        Points[i].second = rk; // Modify the Y-coordinate for fake shares
    }

    // 3
    vector<vector<long long>> KM(s, vector<long long>(t));
    M.clear();
    for (int i = 0; i < s; ++i)
        for (int j = 0; j < t; ++j)
            KM[i][j] = get_rand(MOD);

    for (int i = 0; i < fake; ++i) { // fake DHF value
        vector<long long> r(s + 1);
        r[0] = Points[i].first;
        for (int j = 1; j <= s; ++j) // Only randomize the s auxiliary evaluations
            r[j] = get_rand(MOD);
        M.push_back(r);
    }

    for (int i = fake; i < real + fake; ++i) { // real DHF value
        vector<long long> r(s + 1);
        r[0] = Points[i].first;
        for (int j = 0; j < s; ++j) {
            long long temp = KM[j][0];
            for (int l = 1; l < t; ++l)
                temp = (temp + KM[j][l] * power(r[0], l)) % MOD;
            r[j + 1] = temp;
        }
        M.push_back(r);
    }

    // 4
    for (int i = 0; i < real + fake; ++i) {
        int k = get_rand(fake + real);

        auto T1 = Points[i];
        Points[i] = Points[k];
        Points[k] = T1;

        auto T2 = M[i];
        M[i] = M[k];
        M[k] = T2;
    }
}

inline long long Detection_Algorithm(long long s, long long t, vector<pair<long long, long long>>& Points, vector<vector<long long>>& v) {
    long long m = v.size();
    vector<vector<long long>> M(s + t, vector<long long>(m));

    // 1. Construct the top Vandermonde basis rows
    for (int i = 0; i < t; ++i)
        for (int j = 0; j < m; ++j)
            M[i][j] = power(v[j][0], i); // Assuming power applies % MOD internally

    // 2. Append the s auxiliary polynomial rows
    for (int i = t; i < t + s; ++i)
        for (int j = 0; j < m; ++j)
            M[i][j] = v[j][i - t + 1];

    // 3. Compute Null Space
    vector<vector<long long>> NL = Null_Space(M, false);
    if (NL.empty())
        return -1; // Trivial null space; not enough valid shares

    vector<long long> b;
    vector<vector<long long>> A;

    // 4. Identify valid shares and collect exactly t+1 points
    for (int i = 0; i < m && b.size() < t + 1; ++i) { // FIX: Changed <= to <

        // FIX: Check the entire row in the null space basis, not just the first column
        bool is_valid = false;
        for (int k = 0; k < NL[i].size(); ++k) {
            if (NL[i][k] != 0) {
                is_valid = true;
                break;
            }
        }

        if (is_valid) {
            b.push_back(Points[i].second);
            vector<long long> tv;
            for (int j = 0; j < t + 1; ++j)
                tv.push_back(power(Points[i].first, j));
            A.push_back(tv);
        }
    }

    // Security check: If we somehow didn't find enough valid shares 
    // despite a non-trivial null space (shouldn't happen in theory, but safe to check)
    if (b.size() < t + 1)
        return -1;

    // 5. Reconstruct K
    vector<long long> x = Ax_b(A, b);
    return x[0]; // The secret K is the y-intercept (coefficient of x^0)
}


inline void testing() {
    double avt = 0;
    for (int trial = 1; trial < 1000000000; ++trial) {
        vector<pair<long long, long long>> v;
        vector<vector<long long>> M;
        long long K = rand() % MOD, s = 1001, t = 10, real = 11, fake = 1000; //only when real > t, it is possible to reconstruct K 
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
