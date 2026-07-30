#include <fstream>
#include <time.h>
#include <chrono>
#include <random>

#include <cuda_runtime.h>

#include <tbb/spin_mutex.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

#include <set>
#include <map>
#include <math.h>
#include <atomic>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>

using namespace std;
long long MOD = 10e9; //should be inside of INT32 range. It's GF so MOD must be a prime.
long long primitive;
vector<int> int_inverse;
vector<int> seeds; // primitive ^ seeds[i] = i
vector<int> cubic_z;
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
    ifstream inFile(filename, ios::binary);
    if (!inFile) {
        cerr << "Error opening file for reading." << endl;
        return;
    }
    auto readVector = [&inFile](auto& vec) {
        size_t size;
        inFile.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(decltype(vec[0])));
        };

    inFile.read(reinterpret_cast<char*>(&primitive), sizeof(primitive));

    // Read vectors
    readVector(int_inverse);
    readVector(seeds);
    readVector(MOD_decompose);
    readVector(MOD_divisors);
    readVector(cubic_z); // Read the new cubic table

    size_t outerSize;
    inFile.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));
    ones_roots.resize(outerSize);
    for (auto& vec : ones_roots)
        readVector(vec);

    inFile.close();
    if (!inFile.good()) {
        cerr << "Error occurred during file read." << endl;
        return;
    }
}
void writeData(const string& filename) {
    ofstream outFile(filename, ios::binary | ios::out);
    if (!outFile) {
        cerr << "Error opening file for writing." << endl;
        return;
    }
    outFile.write(reinterpret_cast<char*>(&primitive), sizeof(primitive));

    auto writeVector = [&outFile](const auto& vec) {
        size_t size = vec.size();
        outFile.write(reinterpret_cast<const char*>(&size), sizeof(size));
        outFile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(decltype(vec[0])));
        };

    // Write vectors to file
    writeVector(int_inverse);
    writeVector(seeds);
    writeVector(MOD_decompose);
    writeVector(MOD_divisors);
    writeVector(cubic_z); // Write the new cubic table

    size_t outerSize = ones_roots.size();
    outFile.write(reinterpret_cast<const char*>(&outerSize), sizeof(outerSize));
    for (const auto& vec : ones_roots)
        writeVector(vec);
    outFile.close();
    if (!outFile.good()) {
        cerr << "Error occurred during file write." << endl;
        return;
    }
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
        for (long long i = 1, t = primitive; i < MOD; ++i, t = (t * primitive) % MOD)
            seeds[t] = (int)i;  // primitive ^ seeds[i] = i

        // Create an array of spin mutexes.
        std::unique_ptr<tbb::spin_mutex[]> locks(new tbb::spin_mutex[MOD_divisors.size()]);
        tbb::parallel_for(1LL, MOD, [&](long long i) {
            for (int j = 0; j < (int)MOD_divisors.size() - 1; ++j) {
                if (power(i, MOD_divisors[j]) == 1) {
                    tbb::spin_mutex::scoped_lock lock(locks[j]);
                    ones_roots[MOD_divisors[j]].push_back((int)i);
                }
            }
            });

        // ==========================================
        // Precompute 1-Parameter Cubic Table (All MODs)
        // M z^3 + z + 1 = 0
        // ==========================================
        cubic_z.resize(MOD, 0);
        for (long long z = 1; z < MOD; ++z) {
            long long z3 = (z * z % MOD) * z % MOD;

            // Fast inverse via primitive seeds (O(1)) instead of power()
            long long z3_inv = power(primitive, MOD - 1 - seeds[z3]);

            long long M = (MOD - z - 1) % MOD * z3_inv % MOD;
            if (M < 0) M += MOD;

            cubic_z[M] = (int)z;
        }

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

__global__ void matrixMulModKernel(const long long* A, const long long* B, long long* C, int rowsA, int colsA, int colsB, long long mod) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < rowsA && col < colsB) {
        unsigned long long sum = 0;
        for (int k = 0; k < colsA; ++k) {
            sum += (unsigned long long)A[row * colsA + k] * B[k * colsB + col];
            if (k & 1) sum %= mod;
        }
        C[row * colsB + col] = sum % mod;
    }
}
inline vector<vector<long long>> operator*(const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty()) return {};
    if (a[0].size() != b.size()) {
        ::printf("Matrix Multiplication Error : Matrix size does not match\n");
        ::exit(1);
    }

    int rowsA = a.size();
    int colsA = a[0].size();
    int colsB = b[0].size();

    if (rowsA <= 10 && colsA <= 10 && colsB <= 10) {
        vector<vector<long long>> R(rowsA, vector<long long>(colsB, 0));
        for (int i = 0; i < rowsA; ++i) {
            for (int k = 0; k < colsA; ++k) {
                long long a_ik = a[i][k];
                for (int j = 0; j < colsB; ++j)
                    R[i][j] = (R[i][j] + a_ik * b[k][j]) % MOD;
            }
        }
        return R;
    }

    vector<long long> flat_A(rowsA * colsA);
    vector<long long> flat_B(colsA * colsB);

    for (int i = 0; i < rowsA; ++i)
        copy(a[i].begin(), a[i].end(), flat_A.begin() + i * colsA);
    for (int i = 0; i < colsA; ++i)
        copy(b[i].begin(), b[i].end(), flat_B.begin() + i * colsB);

    size_t bytesA = rowsA * colsA * sizeof(long long);
    size_t bytesB = colsA * colsB * sizeof(long long);
    size_t bytesC = rowsA * colsB * sizeof(long long);

    long long* d_A, * d_B, * d_C;
    cudaMalloc(&d_A, bytesA);
    cudaMalloc(&d_B, bytesB);
    cudaMalloc(&d_C, bytesC);

    cudaMemcpy(d_A, flat_A.data(), bytesA, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, flat_B.data(), bytesB, cudaMemcpyHostToDevice);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((colsB + 15) / 16, (rowsA + 15) / 16);

    // Pass the rowsA, colsA, colsB parameters explicitly to the kernel!
    matrixMulModKernel << <numBlocks, threadsPerBlock >> > (d_A, d_B, d_C, rowsA, colsA, colsB, MOD);

    vector<long long> flat_C(rowsA * colsB);
    cudaMemcpy(flat_C.data(), d_C, bytesC, cudaMemcpyDeviceToHost);

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    vector<vector<long long>> R(rowsA, vector<long long>(colsB));
    for (int i = 0; i < rowsA; ++i)
        copy(flat_C.begin() + i * colsB, flat_C.begin() + (i + 1) * colsB, R[i].begin());

    return R;
}
__global__ void matrixVectorMulModKernel(const long long* A, const long long* B, long long* C, int rowsA, int colsA, long long mod) {
    // 1D grid/block layout since the output is a 1D vector
    int row = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < rowsA) {
        unsigned long long sum = 0;
        for (int j = 0; j < colsA; ++j) {
            sum += (unsigned long long)A[row * colsA + j] * B[j];
            // Match the modulo optimization from your original code
            if (j & 1) sum %= mod;
        }
        C[row] = sum % mod;
    }
}
inline vector<long long> operator*(const vector<vector<long long>>& a, const vector<long long>& b) {
    if (a.empty() || a[0].empty() || b.empty()) return {};
    if (a[0].size() != b.size()) {
        printf("Matrix Vector Multiplication Error : Matrix and Vector's size do not match\n");
        exit(1);
    }

    int rowsA = a.size();
    int colsA = a[0].size();

    // CPU fallback for small workloads to avoid CUDA launch overhead
    if (rowsA <= 20 && colsA <= 20) {
        vector<long long> R(rowsA, 0);
        for (int i = 0; i < rowsA; ++i) {
            unsigned long long sum = 0;
            for (int j = 0; j < colsA; ++j) {
                sum += (unsigned long long)a[i][j] * b[j];
                if (j & 1) sum %= MOD;
            }
            R[i] = sum % MOD;
        }
        return R;
    }

    // Flatten 2D matrix A to 1D
    vector<long long> flat_A(rowsA * colsA);
    for (int i = 0; i < rowsA; ++i) {
        copy(a[i].begin(), a[i].end(), flat_A.begin() + i * colsA);
    }

    // Memory sizing
    size_t bytesA = rowsA * colsA * sizeof(long long);
    size_t bytesB = colsA * sizeof(long long);
    size_t bytesC = rowsA * sizeof(long long);

    // Allocate Device memory
    long long* d_A, * d_B, * d_C;
    cudaMalloc(&d_A, bytesA);
    cudaMalloc(&d_B, bytesB);
    cudaMalloc(&d_C, bytesC);

    // Copy Host to Device
    cudaMemcpy(d_A, flat_A.data(), bytesA, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, b.data(), bytesB, cudaMemcpyHostToDevice); // b is already flat

    // 1D Thread/Block configuration
    int threads = 256;
    int blocks = (rowsA + threads - 1) / threads;

    // Launch kernel
    matrixVectorMulModKernel << <blocks, threads >> > (d_A, d_B, d_C, rowsA, colsA, MOD);

    // Copy Device to Host
    vector<long long> R(rowsA);
    cudaMemcpy(R.data(), d_C, bytesC, cudaMemcpyDeviceToHost);

    // Free Device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

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
            ::printf("Matrix Chop Error : Block bounds exceed matrix dimensions\n");
            ::exit(1);
        }

        vector<vector<long long>> block(block_size, vector<long long>(block_size));
        for (size_t j = 0; j < block_size; ++j) {
            // Source Start: F[row].begin() + column_offset
            // Source End: F[row].begin() + column_offset + block_size
            // Destination: block[j].begin()
            copy(
                F[j + p].begin() + p,
                F[j + p].begin() + p + block_size,
                block[j].begin()
            );
        }
        M.push_back(move(block));
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

// Kernel to swap two rows in parallel
__global__ void swapRowsKernel(long long* A, int n, int row1, int row2) {
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c < n) {
        long long temp = A[row1 * n + c];
        A[row1 * n + c] = A[row2 * n + c];
        A[row2 * n + c] = temp;
    }
}

// Kernel to multiply a row by the pivot's modular inverse
__global__ void normalizeRowKernel(long long* A, int n, int row, int col, long long inv, long long mod) {
    int c = col + blockIdx.x * blockDim.x + threadIdx.x;
    if (c < n) {
        A[row * n + c] = (A[row * n + c] * inv) % mod;
    }
}

// Kernel to eliminate all rows below the current pivot row
__global__ void eliminateRowsKernel(long long* A, int m, int n, int pivot_row, int col, long long mod) {
    // Grid maps: X -> columns (starting from col), Y -> rows (starting from pivot_row + 1)
    int r = (pivot_row + 1) + blockIdx.y * blockDim.y + threadIdx.y;
    int c = col + blockIdx.x * blockDim.x + threadIdx.x;

    if (r < m && c < n) {
        long long factor = A[r * n + col];
        if (factor != 0) { // Optimization: Skip if already 0
            long long pivot_val = A[pivot_row * n + c];
            long long sub = (pivot_val * factor) % mod;
            long long new_val = A[r * n + c] - sub;

            if (new_val < 0) new_val += mod;
            A[r * n + c] = new_val;
        }
    }
}
inline size_t matrix_rank_gpu(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D[0].empty()) return 0;

    int m = A_2D.size();
    int n = A_2D[0].size();
    size_t bytes = m * n * sizeof(long long);

    // Allocate Unified Memory (accessible by both CPU and GPU)
    long long* d_A = nullptr;
    cudaMallocManaged(&d_A, bytes);

    // Copy the 2D vector data into our Unified Memory buffer
    for (int i = 0; i < m; ++i)
        copy(A_2D[i].begin(), A_2D[i].end(), d_A + i * n);

    size_t rank = 0;
    int row = 0; // Tracks the current pivot row

    for (int col = 0; col < n && row < m; ++col) {
        // 1. Pivot Search (CPU-side on Unified Memory)
        // We must sync the GPU first to ensure previous eliminations are complete
        cudaDeviceSynchronize();

        int pivot_row = row;
        while (pivot_row < m && d_A[pivot_row * n + col] == 0) {
            pivot_row++;
        }

        if (pivot_row == m) continue; // All zeros in this column below 'row'

        // 2. Row Swapping (GPU-side)
        if (pivot_row != row) {
            int threads = 256;
            int blocks = (n + threads - 1) / threads;
            swapRowsKernel << <blocks, threads >> > (d_A, n, row, pivot_row);
        }

        // 3. Row Normalization (GPU-side)
        // Calculate the inverse on the CPU using your existing 'inverse()' function
        long long inv = inverse(d_A[row * n + col]);

        int threads_norm = 256;
        int blocks_norm = ((n - col) + threads_norm - 1) / threads_norm;
        normalizeRowKernel << <blocks_norm, threads_norm >> > (d_A, n, row, col, inv, MOD);

        // 4. Row Elimination (GPU-side)
        int active_rows = m - (row + 1);
        int active_cols = n - col;

        if (active_rows > 0 && active_cols > 0) {
            dim3 threadsPerBlock(16, 16);
            dim3 numBlocks((active_cols + 15) / 16, (active_rows + 15) / 16);

            eliminateRowsKernel << <numBlocks, threadsPerBlock >> > (d_A, m, n, row, col, MOD);
        }

        row++;
        rank++;
    }

    // Wait for final GPU operations to finish before freeing memory
    cudaDeviceSynchronize();
    cudaFree(d_A);

    return rank;
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
inline void matrix_diagonalize_2x2(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D) {
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
inline void matrix_diagonalize_3x3(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D) {
    S.assign(3, vector<long long>(3, 0));
    D.assign(3, vector<long long>(3, 0));

    // Fast mapping for matrix components
    long long A00 = A[0][0], A01 = A[0][1], A02 = A[0][2];
    long long A10 = A[1][0], A11 = A[1][1], A12 = A[1][2];
    long long A20 = A[2][0], A21 = A[2][1], A22 = A[2][2];

    // 1. Trace (c2), Sum of Principal Minors (c1), and Determinant (c0)
    long long c2 = (A00 + A11 + A22) % MOD;

    long long m00 = (A11 * A22 - A12 * A21) % MOD;
    long long m11 = (A00 * A22 - A02 * A20) % MOD;
    long long m22 = (A00 * A11 - A01 * A10) % MOD;

    long long c1 = (m00 + m11 + m22) % MOD;
    if (c1 < 0) c1 += MOD;

    // Sub-minors used for c0
    long long c0 = (A00 * m00 - A01 * ((A10 * A22 - A12 * A20) % MOD) + A02 * ((A10 * A21 - A11 * A20) % MOD)) % MOD;
    if (c0 < 0) c0 += MOD;

    long long lambdas[3];
    long long inv2 = inverse(2);

    // 2. Branch depending on whether primitive 3rd roots of unity exist in F_p
    if ((MOD - 1) % 3 == 0) {
        long long inv3 = inverse(3);
        long long inv27 = inverse(27);

        long long P = (c1 - c2 * c2 % MOD * inv3) % MOD;
        if (P < 0) P += MOD;

        long long c2_3 = c2 * c2 % MOD * c2 % MOD;
        long long Q = (-c0 + c1 * c2 % MOD * inv3 - 2 * c2_3 % MOD * inv27) % MOD;
        if (Q < 0) Q += MOD;

        long long P3 = P * P % MOD * P % MOD;
        long long Delta_quad = (Q * Q % MOD + 4 * P3 % MOD * inv27) % MOD;
        if (Delta_quad < 0) Delta_quad += MOD;

        long long S_res = Delta_quad != 0 ? power(primitive, seeds[Delta_quad] >> 1) : 0;

        long long t1 = (MOD - Q + S_res) % MOD * inv2 % MOD;
        long long t2 = (MOD * 2 - Q - S_res) % MOD * inv2 % MOD;

        long long u = t1 != 0 ? power(primitive, seeds[t1] / 3) : 0;
        long long v = 0;

        if (P == 0) {
            if (t2 != 0) v = power(primitive, seeds[t2] / 3);
        }
        else {
            v = (MOD - P) * inv3 % MOD * inverse(u) % MOD;
        }

        long long omega = ones_roots[3][1];
        long long omega2 = ones_roots[3][2];
        long long shift = c2 * inv3 % MOD;

        lambdas[0] = (u + v + shift) % MOD;
        lambdas[1] = (u * omega % MOD + v * omega2 % MOD + shift) % MOD;
        lambdas[2] = (u * omega2 % MOD + v * omega % MOD + shift) % MOD;
    }
    else {
        long long lambda1 = -1;
        auto check_root = [&](long long r) {
            long long r2 = r * r % MOD;
            return (r2 * r - c2 * r2 + c1 * r - c0) % MOD == 0;
            };

        if (check_root(0)) {
            lambda1 = 0;
        }
        else {
            long long a = 1, half_p = (MOD - 1) >> 1;
            while (true) {
                long long res_u = 0, res_v = 0, res_w = 1;
                long long base_u = 0, base_v = 1, base_w = a;
                long long exp = half_p;

                while (exp > 0) {
                    if (exp & 1) {
                        long long C4 = res_u * base_u % MOD;
                        long long C3 = (res_u * base_v + res_v * base_u) % MOD;
                        long long C2 = (res_u * base_w + res_v * base_v + res_w * base_u) % MOD;
                        long long C1 = (res_v * base_w + res_w * base_v) % MOD;
                        long long C0 = res_w * base_w % MOD;

                        C3 = (C3 + C4 * c2) % MOD; C2 = (C2 - C4 * c1) % MOD; C1 = (C1 + C4 * c0) % MOD;

                        res_u = (C2 + C3 * c2) % MOD; if (res_u < 0) res_u += MOD;
                        res_v = (C1 - C3 * c1) % MOD; if (res_v < 0) res_v += MOD;
                        res_w = (C0 + C3 * c0) % MOD; if (res_w < 0) res_w += MOD;
                    }
                    long long C4 = base_u * base_u % MOD;
                    long long C3 = 2 * base_u * base_v % MOD;
                    long long C2 = (2 * base_u * base_w + base_v * base_v) % MOD;
                    long long C1 = 2 * base_v * base_w % MOD;
                    long long C0 = base_w * base_w % MOD;

                    C3 = (C3 + C4 * c2) % MOD; C2 = (C2 - C4 * c1) % MOD; C1 = (C1 + C4 * c0) % MOD;

                    base_u = (C2 + C3 * c2) % MOD; if (base_u < 0) base_u += MOD;
                    base_v = (C1 - C3 * c1) % MOD; if (base_v < 0) base_v += MOD;
                    base_w = (C0 + C3 * c0) % MOD; if (base_w < 0) base_w += MOD;
                    exp >>= 1;
                }

                res_w = (res_w + MOD - 1) % MOD;
                if (res_u == 0 && res_v == 0) { a++; continue; }

                if (res_u == 0) {
                    long long cand = (MOD - res_w) * inverse(res_v) % MOD;
                    if (check_root(cand)) { lambda1 = cand; break; }
                }
                else {
                    long long inv_u = inverse(res_u);
                    long long V = res_v * inv_u % MOD;
                    long long W = res_w * inv_u % MOD;

                    long long A_div = (MOD * 2 - c2 - V) % MOD;
                    long long R1 = (c1 - W - A_div * V) % MOD; if (R1 < 0) R1 += MOD;
                    long long R0 = (-c0 - A_div * W) % MOD;    if (R0 < 0) R0 += MOD;

                    if (R1 == 0 && R0 == 0) {
                        long long cand = (V + c2) % MOD;
                        if (check_root(cand)) { lambda1 = cand; break; }
                    }
                    else if (R1 != 0) {
                        long long cand = (MOD - R0) * inverse(R1) % MOD;
                        if (check_root(cand)) { lambda1 = cand; break; }
                    }
                }
                a++;
            }
        }

        // Deflate
        long long A_coef = (lambda1 - c2 + MOD) % MOD;
        long long B_coef = (c1 + lambda1 * A_coef) % MOD; if (B_coef < 0) B_coef += MOD;

        long long discriminant = (A_coef * A_coef - 4 * B_coef) % MOD;
        if (discriminant < 0) discriminant += MOD;

        long long sq = discriminant != 0 ? power(primitive, seeds[discriminant] >> 1) : 0;

        lambdas[0] = lambda1;
        lambdas[1] = (MOD - A_coef + sq) % MOD * inv2 % MOD;
        lambdas[2] = (MOD - A_coef - sq + MOD) % MOD * inv2 % MOD;
    }

    D[0][0] = lambdas[0]; D[1][1] = lambdas[1]; D[2][2] = lambdas[2];

    // 3. Construct Eigenvectors Fast Column Mapping
    bool used[3] = { false, false, false };
    for (int i = 0; i < 3; i++) {
        if (used[i]) continue;
        long long lambda = lambdas[i];

        int indices[3], m = 0;
        for (int j = i; j < 3; j++) {
            if (lambdas[j] == lambda) { indices[m++] = j; used[j] = true; }
        }

        // Direct variables (No 3x3 allocations per iteration)
        long long M00 = A00 - lambda, M01 = A01, M02 = A02;
        long long M10 = A10, M11 = A11 - lambda, M12 = A12;
        long long M20 = A20, M21 = A21, M22 = A22 - lambda;
        M00 = M00 < 0 ? M00 + MOD : M00;
        M11 = M11 < 0 ? M11 + MOD : M11;
        M22 = M22 < 0 ? M22 + MOD : M22;

        if (m == 1) {
            long long v0 = (M01 * M12 - M02 * M11) % MOD;
            long long v1 = (M02 * M10 - M00 * M12) % MOD;
            long long v2 = (M00 * M11 - M01 * M10) % MOD;

            if (v0 == 0 && v1 == 0 && v2 == 0) {
                v0 = (M11 * M22 - M12 * M21) % MOD;
                v1 = (M12 * M20 - M10 * M22) % MOD;
                v2 = (M10 * M21 - M11 * M20) % MOD;
            }
            if (v0 == 0 && v1 == 0 && v2 == 0) {
                v0 = (M01 * M22 - M02 * M21) % MOD;
                v1 = (M02 * M20 - M00 * M22) % MOD;
                v2 = (M00 * M21 - M01 * M20) % MOD;
            }

            S[0][indices[0]] = v0 < 0 ? v0 + MOD : v0;
            S[1][indices[0]] = v1 < 0 ? v1 + MOD : v1;
            S[2][indices[0]] = v2 < 0 ? v2 + MOD : v2;
        }
        else if (m == 2) {
            long long r0 = 0, r1 = 0, r2 = 0;
            if (M00 || M01 || M02) { r0 = M00; r1 = M01; r2 = M02; }
            else if (M10 || M11 || M12) { r0 = M10; r1 = M11; r2 = M12; }
            else { r0 = M20; r1 = M21; r2 = M22; }

            if (r0 != 0) {
                S[0][indices[0]] = (MOD - r1) % MOD; S[1][indices[0]] = r0; S[2][indices[0]] = 0;
                S[0][indices[1]] = (MOD - r2) % MOD; S[1][indices[1]] = 0;  S[2][indices[1]] = r0;
            }
            else if (r1 != 0) {
                S[0][indices[0]] = 1; S[1][indices[0]] = 0; S[2][indices[0]] = 0;
                S[0][indices[1]] = 0; S[1][indices[1]] = (MOD - r2) % MOD; S[2][indices[1]] = r1;
            }
            else {
                S[0][indices[0]] = 1; S[1][indices[0]] = 0; S[2][indices[0]] = 0;
                S[0][indices[1]] = 0; S[1][indices[1]] = 1; S[2][indices[1]] = 0;
            }
        }
        else if (m == 3) {
            S[0][indices[0]] = 1; S[1][indices[0]] = 0; S[2][indices[0]] = 0;
            S[0][indices[1]] = 0; S[1][indices[1]] = 1; S[2][indices[1]] = 0;
            S[0][indices[2]] = 0; S[1][indices[2]] = 0; S[2][indices[2]] = 1;
        }
    }
}
inline void matrix_diagonalize_4x4(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D) {
    S.assign(4, vector<long long>(4, 0));
    D.assign(4, vector<long long>(4, 0));

    // Fast mapping for matrix components
    long long A00 = A[0][0], A01 = A[0][1], A02 = A[0][2], A03 = A[0][3];
    long long A10 = A[1][0], A11 = A[1][1], A12 = A[1][2], A13 = A[1][3];
    long long A20 = A[2][0], A21 = A[2][1], A22 = A[2][2], A23 = A[2][3];
    long long A30 = A[3][0], A31 = A[3][1], A32 = A[3][2], A33 = A[3][3];

    // Precompute constant inverses
    long long inv2 = inverse(2);
    long long inv3 = inverse(3);
    long long inv4 = inverse(4);
    long long inv8 = inverse(8);
	long long inv16 = inverse(16);
    long long inv27 = inverse(27);
    long long inv256 = inverse(256%MOD);

    // 1. Calculate Characteristic Polynomial Coefficients (x^4 + a*x^3 + b*x^2 + c*x + d = 0)
    long long S1 = (A00 + A11 + A22 + A33) % MOD;

    auto minor2 = [&](long long m11, long long m12, long long m21, long long m22) {
        long long val = (m11 * m22 - m12 * m21) % MOD;
        return val < 0 ? val + MOD : val;
        };
    long long S2 = (minor2(A00, A01, A10, A11) + minor2(A00, A02, A20, A22) + minor2(A00, A03, A30, A33) +
        minor2(A11, A12, A21, A22) + minor2(A11, A13, A31, A33) + minor2(A22, A23, A32, A33)) % MOD;

    auto minor3 = [&](long long m11, long long m12, long long m13,
        long long m21, long long m22, long long m23,
        long long m31, long long m32, long long m33) {
            long long val = m11 * minor2(m22, m23, m32, m33) - m12 * minor2(m21, m23, m31, m33) + m13 * minor2(m21, m22, m31, m32);
            val %= MOD;
            return val < 0 ? val + MOD : val;
        };
    long long S3 = (minor3(A00, A01, A02, A10, A11, A12, A20, A21, A22) +
        minor3(A00, A01, A03, A10, A11, A13, A30, A31, A33) +
        minor3(A00, A02, A03, A20, A22, A23, A30, A32, A33) +
        minor3(A11, A12, A13, A21, A22, A23, A31, A32, A33)) % MOD;

    long long S4 = (A00 * minor3(A11, A12, A13, A21, A22, A23, A31, A32, A33) -
        A01 * minor3(A10, A12, A13, A20, A22, A23, A30, A32, A33) +
        A02 * minor3(A10, A11, A13, A20, A21, A23, A30, A31, A33) -
        A03 * minor3(A10, A11, A12, A20, A21, A22, A30, A31, A32)) % MOD;
    if (S4 < 0) S4 += MOD;

    long long a = MOD - S1;
    long long b = S2;
    long long c = MOD - S3;
    long long d = S4;

    // 2. Depress the Quartic (y^4 + p*y^2 + q*y + r = 0) where x = y - a/4
    long long a2 = (a * a) % MOD;
    long long a3 = (a2 * a) % MOD;
    long long a4 = (a3 * a) % MOD;

    long long p = (b - 3 * a2 % MOD * inv8) % MOD;
    if (p < 0) p += MOD;

    long long q = (c - a * b % MOD * inv2 % MOD + a3 * inv8) % MOD;
    if (q < 0) q += MOD;

    long long r_val = (d - a * c % MOD * inv4 % MOD + a2 * b % MOD * inv16 % MOD - 3 * a4 % MOD * inv256 % MOD) % MOD;
    if (r_val < 0) r_val += MOD;

    // 3. Resolvent Cubic (z^3 + Ac*z^2 + Bc*z + Cc = 0)
    long long p2 = (p * p) % MOD;
    long long q2 = (q * q) % MOD;

    long long Ac = p;
    long long Bc = (p2 * inv4 - r_val) % MOD;
    if (Bc < 0) Bc += MOD;
    long long Cc = (MOD - q2 * inv8 % MOD) % MOD;

    // Depress the Cubic (t^3 + Pc*t + Qc = 0) where z = t - Ac/3
    long long Pc = (Bc - Ac * Ac % MOD * inv3) % MOD;
    if (Pc < 0) Pc += MOD;

    long long Qc = (Cc + 2 * Ac * Ac % MOD * Ac % MOD * inv27 % MOD - Ac * Bc % MOD * inv3) % MOD;
    if (Qc < 0) Qc += MOD;

    // 4. Solve the Cubic in True O(1)
    long long t = 0;
    if (Pc == 0) { // t^3 = -Qc
        long long K = MOD - Qc;
        if (K != 0) {
            long long s = seeds[K];
            if ((MOD - 1) % 3 == 0) t = power(primitive, s / 3);
            else {
                long long inv3_phi = (2 * MOD - 1) / 3;
                t = power(primitive, (s * inv3_phi) % (MOD - 1));
            }
        }
    }
    else if (Qc == 0) {
        t = 0;
    }
    else { // 1-Parameter Table Lookup
        long long P3 = Pc * Pc % MOD * Pc % MOD;
        long long P3_inv = inverse(P3);
        long long M = Qc * Qc % MOD * P3_inv % MOD;
        long long z_val = cubic_z[M];
        t = z_val * Qc % MOD * inverse(Pc) % MOD;
    }

    long long z = (t - Ac * inv3) % MOD;
    if (z < 0) z += MOD;

    // 5. Split Quartic into Two Quadratics
    long long lambdas[4];
    long long z2_val = (2 * z) % MOD;

    if (z2_val == 0) { // Bi-quadratic case (q = 0)
        long long delta_u = (p * p - 4 * r_val) % MOD;
        if (delta_u < 0) delta_u += MOD;

        long long sq_u = delta_u != 0 ? power(primitive, seeds[delta_u] >> 1) : 0;
        long long u1 = (MOD - p + sq_u) % MOD * inv2 % MOD;
        long long u2 = (MOD - p - sq_u + MOD) % MOD * inv2 % MOD;

        long long y1 = 0, y2 = 0, y3 = 0, y4 = 0;
        if (u1 != 0) { long long s1 = power(primitive, seeds[u1] >> 1); y1 = s1; y2 = MOD - s1; }
        if (u2 != 0) { long long s2 = power(primitive, seeds[u2] >> 1); y3 = s2; y4 = MOD - s2; }

        lambdas[0] = y1; lambdas[1] = y2; lambdas[2] = y3; lambdas[3] = y4;
    }
    else { // Standard Ferrari split
        long long W = power(primitive, seeds[z2_val] >> 1);
        long long inv2W = inverse((2 * W) % MOD);

        long long term1 = (p * inv2 + z) % MOD;
        long long term2 = (q * inv2W) % MOD;

        long long C1 = term1 - term2; if (C1 < 0) C1 += MOD;
        long long C2 = term1 + term2; if (C2 < 0) C2 += MOD;

        // Quadratic 1: y^2 + Wy + C1 = 0
        long long delta1 = (W * W - 4 * C1) % MOD; if (delta1 < 0) delta1 += MOD;
        long long sq1 = delta1 != 0 ? power(primitive, seeds[delta1] >> 1) : 0;
        lambdas[0] = (MOD - W + sq1) % MOD * inv2 % MOD;
        lambdas[1] = (MOD - W - sq1 + MOD) % MOD * inv2 % MOD;

        // Quadratic 2: y^2 - Wy + C2 = 0
        long long delta2 = (W * W - 4 * C2) % MOD; if (delta2 < 0) delta2 += MOD;
        long long sq2 = delta2 != 0 ? power(primitive, seeds[delta2] >> 1) : 0;
        lambdas[2] = (W + sq2) % MOD * inv2 % MOD;
        lambdas[3] = (W - sq2 + MOD) % MOD * inv2 % MOD;
    }

    // Shift Eigenvalues back: lambda = y - a/4
    long long shift = (MOD - a * inv4 % MOD) % MOD;
    for (int i = 0; i < 4; i++) {
        lambdas[i] = (lambdas[i] + shift) % MOD;
        D[i][i] = lambdas[i];
    }

    // 6. Eigenspace (Null Space) Construction via Fixed-Size Gauss-Jordan
    bool used[4] = { false, false, false, false };
    for (int i = 0; i < 4; i++) {
        if (used[i]) continue;
        long long lambda = lambdas[i];

        int indices[4], m = 0;
        for (int j = i; j < 4; j++)
            if (lambdas[j] == lambda) { indices[m++] = j; used[j] = true; }

        long long M[4][4];
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) M[r][c] = A[r][c];
            M[r][r] = (M[r][r] - lambda + MOD) % MOD;
        }

        int pivot_col[4] = { -1, -1, -1, -1 };
        int row = 0;
        for (int col = 0; col < 4 && row < 4; ++col) {
            int pivot = row;
            while (pivot < 4 && M[pivot][col] == 0) pivot++;
            if (pivot == 4) continue;

            if (pivot != row)
                for (int c = col; c < 4; ++c) swap(M[row][c], M[pivot][c]);

            long long inv_val = inverse(M[row][col]);
            for (int c = col; c < 4; ++c) M[row][c] = (M[row][c] * inv_val) % MOD;

            for (int r = 0; r < 4; ++r) {
                if (r != row && M[r][col] != 0) {
                    long long factor = M[r][col];
                    for (int c = col; c < 4; ++c) {
                        M[r][c] = (M[r][c] - factor * M[row][c]) % MOD;
                        if (M[r][c] < 0) M[r][c] += MOD;
                    }
                }
            }
            pivot_col[row] = col;
            row++;
        }

        bool is_free[4] = { true, true, true, true };
        for (int r = 0; r < row; ++r) is_free[pivot_col[r]] = false;

        int m_count = 0;
        for (int col = 0; col < 4; ++col) {
            if (is_free[col] && m_count < m) {
                S[col][indices[m_count]] = 1;
                for (int r = 0; r < row; ++r)
                    S[pivot_col[r]][indices[m_count]] = (MOD - M[r][col]) % MOD;
                m_count++;
            }
        }

        if (m_count < m) {
            printf("Matrix Diagonalize Error : 4x4 matrix is defective (not diagonalizable)\n\n");
            exit(1);
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
void matrix_diagonalize_krylov(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, long long MOD);
inline void matrix_diagonalize_henry(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    int n = (int)A.size();
    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    vector<vector<long long>> AP_1 = matrix_power(A, MOD - 1);
    vector<vector<long long>> ZN1 = Null_Space(AP_1 - I_n(n), Orth);
    int eigvec_count = ZN1.empty() ? 0 : (int)ZN1[0].size();
    if (eigvec_count > 0)
        for (int j = 0; j < n; ++j)
            copy(ZN1[j].begin(), ZN1[j].end(), S[j].begin());
    vector<vector<long long>> ZN2 = Null_Space(AP_1, Orth);
    if (!ZN2.empty())
        for (int j = 0; j < n; ++j)
            copy(ZN2[j].begin(), ZN2[j].end(), S[j].begin() + eigvec_count);

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

                    int N = M[m_idx].size();
                    if (N >= 1 && N <= 4) {
                        vector<vector<long long>> D, S;
                        if (N == 1) {
                            D = { { M[m_idx][0][0] } };
                            S = { { 1 } };
                        }
                        else if (N == 2)
                            matrix_diagonalize_2x2(M[m_idx], S, D);
                        else if (N == 3)
                            matrix_diagonalize_3x3(M[m_idx], S, D);
                        else if (N == 4)
                            matrix_diagonalize_4x4(M[m_idx], S, D);
                        for (int r = 0; r < N; ++r) {
                            local_result.new_matrices.push_back({ {D[r][r]} });
                            local_result.new_FEs.push_back(D[r][r]);
                            for (int c = 0; c < N; ++c)
                                ST[local_stp + r][local_stp + c] = S[r][c];
                        }
                        continue;
                    }

                    int m_size = M[m_idx].size();
                    int roots_size = ones_roots[MOD_decompose[pi]].size();
                    if (roots_size > log(m_size)) {
                    //if (roots_size > m_size) {
                            //printf("big(%d %d)",roots_size,m_size);
                        vector<vector<long long>> S_krylov, D_krylov;
                        matrix_diagonalize_krylov(M[m_idx], S_krylov, D_krylov, MOD);
                        for (int i = 0; i < m_size; ++i) {
                            local_result.new_matrices.push_back({ {D_krylov[i][i]} });
                            local_result.new_FEs.push_back(D_krylov[i][i]);
                            for (int r = 0; r < m_size; ++r)
                                ST[local_stp + r][local_stp + i] = S_krylov[r][i];
                        }
                        continue; // Skip the rest of the loop because it is already fully diagonalized
                    }

                    vector<vector<long long>> PM = matrix_power(M[m_idx], powC);
                    vector<vector<long long>> St(m_size, vector<long long>(m_size, 0));
                    vector<int> eigspace_dim;
                    int current_col = 0, local_eigvec_count = 0;

                    long long seed = seeds[FE[m_idx]] * inverse(MOD_decompose[pi]) % MOD;
                    long long seed2 = power(primitive, seed);

                    struct RootResult {
                        long long candidate;
                        vector<vector<long long>> ZN;
                    };
                    vector<RootResult> valid_roots;






                    //atomic<int> local_eigvec_count{ 0 };

                    //tbb::parallel_for(tbb::blocked_range<int>(0, roots_size),
                    //    [&](const tbb::blocked_range<int>& inner_r) {
                    //        for (int i = inner_r.begin(); i != inner_r.end(); ++i) {
                    //            if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;

                    //            long long candidate = seed2 * ones_roots[MOD_decompose[pi]][i] % MOD;
                    //            vector<vector<long long>> query = PM;
                    //            for (int j = 0; j < query.size(); ++j) {
                    //                query[j][j] -= candidate;
                    //                if (query[j][j] < 0) query[j][j] += MOD;
                    //            }

                    //            vector<vector<long long>> ZN = Null_Space(query, Orth);

                    //            if (!ZN.empty()) {
                    //                valid_roots.push_back({ candidate, ZN });
                    //                local_eigvec_count.fetch_add(ZN[0].size(), std::memory_order_relaxed);
                    //            }
                    //        }
                    //    });

                    





                    // Setup all candidates
                    vector<long long> candidates(roots_size);
                    for (int i = 0; i < roots_size; ++i)
                        candidates[i] = seed2 * ones_roots[MOD_decompose[pi]][i] % MOD;

                    vector<long long> found_eigenvalues; // Track found roots globally

                    // Helper to compute: v_out = (PM * v_in) - (c * v_in)
                    auto apply_PM_minus_cI = [&](const vector<long long>& v_in, long long c) -> pair<vector<long long>, bool> {
                        vector<long long> v_out = PM * v_in;
                        bool is_zero = true;
                        for (int r = 0; r < m_size; ++r) {
                            long long val = (v_out[r] - c * v_in[r]) % MOD;
                            if (val < 0) val += MOD;
                            v_out[r] = val;
                            if (val != 0) is_zero = false;
                        }
                        return { move(v_out), is_zero };
                        };

                    // Helper to extract eigenspace and update global state
                    auto extract_and_add_eigenvalue = [&](long long found_val) {
                        vector<vector<long long>> ZN = Null_Space(PM - I_n(PM.size(), found_val), Orth);
                        if (!ZN.empty()) {
                            valid_roots.push_back({ found_val, ZN });
                            local_eigvec_count += ZN[0].size();
                            found_eigenvalues.push_back(found_val);
                        }
                        };

                    // Recursive Divide and Conquer for Eigenvalues
                    auto find_eigenvalues_rec = [&](auto& self, vector<long long> v, int start, int end) -> void {
                        if (local_eigvec_count >= m_size) return;

                        // Base case: if vector is zero
                        bool v_is_zero = true;
                        for (long long x : v) if (x != 0) { v_is_zero = false; break; }
                        if (v_is_zero) return;

                        if (start >= end) return;
                        if (start == end - 1) { // Base case: 1 candidate left
                            extract_and_add_eigenvalue(candidates[start]);
                            return;
                        }
						printf("%d ", end - start);

                        int mid = start + ((end - start) >> 1);

                        // 1. Apply L (start to mid) to isolate R-components
                        vector<long long> v_R = v;
                        for (int i = start; i < mid; ++i) {
                            pair<vector<long long>, bool> res = apply_PM_minus_cI(v_R, candidates[i]);
                            if (res.second) {
                                // MASSIVE OPTIMIZATION: Vector annihilated! 
                                // R is completely empty, and the rest of L (i+1 to mid) is empty.
                                extract_and_add_eigenvalue(candidates[i]);
                                pair<vector<long long>, bool> deflated_v = apply_PM_minus_cI(v, candidates[i]);
                                self(self, deflated_v.first, start, i);
                                return;
                            }
                            v_R = move(res.first);
                        }

                        // 2. Search R (mid to end) using the isolated vector
                        int count_before = found_eigenvalues.size();
                        self(self, v_R, mid, end);
                        int count_after = found_eigenvalues.size();

                        // 3. Deflate newly found eigenvalues from R
                        vector<long long> v_L = move(v);
                        for (int idx = count_before; idx < count_after; ++idx) {
                            pair<vector<long long>, bool> res = apply_PM_minus_cI(v_L, found_eigenvalues[idx]);
                            v_L = move(res.first);
                        }

                        // 4. Search L (start to mid)
                        self(self, v_L, start, mid);
                        };

                    // Iterate through standard basis vectors
                    for (int basis_idx = 0; basis_idx < m_size && local_eigvec_count < m_size; ++basis_idx) {
                        vector<long long> v(m_size, 0);
                        v[basis_idx] = 1;

                        // PRE-CHECK: Is v already a pure eigenvector? (Collinearity check: PM * v == c * v)
                        vector<long long> PM_v = PM * v;
                        long long c = -1;
                        bool is_eigenvector = true;

                        for (int i = 0; i < m_size; ++i) {
                            if (v[i] != 0) {
                                if (c == -1) // Find the scalar c from the first non-zero element
                                    c = (PM_v[i] * inverse(v[i])) % MOD;
                                long long expected = (c * v[i]) % MOD;
                                if (PM_v[i] != expected) {
                                    is_eigenvector = false;
                                    break;
                                }
                            }
                            else if (PM_v[i] != 0) {
                                is_eigenvector = false;
                                break;
                            }
                        }

                        if (is_eigenvector) {
                            // v is a pure eigenvector! Extract it and skip the heavy recursion entirely.
                            extract_and_add_eigenvalue(c);
                            continue;
                        }

                        // Not a pure eigenvector, run the recursive divide-and-conquer
                        find_eigenvalues_rec(find_eigenvalues_rec, v, 0, roots_size);
                        break;
                    }








                    for (auto& root : valid_roots) {
                        if (current_col >= m_size) break;

                        local_result.new_FEs.push_back(root.candidate);
                        eigspace_dim.push_back((int)root.ZN[0].size());

                        for (int j = 0; j < root.ZN.size(); ++j)
                            copy(root.ZN[j].begin(), root.ZN[j].end(), St[j].begin() + current_col);
                        current_col += root.ZN[0].size();
                    }
					//printf("%d", (int)eigspace_dim.size());

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
            D[Di][Di] = M[mat_i][i][i];
    vector<vector<long long>> St_final = I_n((int)S.size());
    for (int i = 0; i < Ss.size(); ++i)
        copy(Ss[i].begin(), Ss[i].end(), St_final[i].begin());

    S = S * St_final;
}





















inline vector<vector<long long>> generate_test_matrix(int N, long long MOD) {
    // 1. Thread-Safe RNG: Pre-generate all random numbers sequentially.
    // Memory allocation and a simple 2D sequential loop is incredibly fast 
    // and completely sidesteps any RNG race conditions.
    vector<vector<long long>> R(N, vector<long long>(N));
    vector<long long> diag(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            R[i][j] = get_rand(MOD);
        }
        diag[i] = get_rand(MOD); // These will be our eigenvalues
    }

    vector<vector<long long>> L(N, vector<long long>(N, 0));
    vector<vector<long long>> U(N, vector<long long>(N, 0));

    // 2. Fully Parallel Construction (No sequential outer loops!)
    tbb::parallel_for(tbb::blocked_range<int>(0, N), [&](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i != r.end(); ++i) {
            L[i][i] = 1; // Unit lower diagonal guarantees det(L) = 1

            for (int j = 0; j < i; ++j) L[i][j] = R[i][j];
            for (int j = i; j < N; ++j) {
                U[i][j] = R[i][j];
                // Guarantee det(U) != 0. (If random hits 0, safely default to 1)
                if (i == j && U[i][i] == 0) U[i][i] = 1;
            }
        }
        });

    // 3. Fully Parallel Matrix Multiplication: S = L * U
    // Optimized: Because L and U are triangular, we only need to loop 'k' up to min(i, j)
    vector<vector<long long>> S(N, vector<long long>(N, 0));
    tbb::parallel_for(tbb::blocked_range<int>(0, N), [&](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i != r.end(); ++i) {
            for (int j = 0; j < N; ++j) {
                unsigned long long local_sum = 0;
                int limit = std::min(i, j);
                for (int k = 0; k <= limit; ++k) {
                    local_sum += (unsigned long long)L[i][k] * U[k][j];
                    if (local_sum >> 61) local_sum %= MOD; // Deferred modulo trick
                }
                S[i][j] = local_sum % MOD;
            }
        }
        });

    // 4. Fully Parallel Column Scaling: SD = S * D
    vector<vector<long long>> SD(N, vector<long long>(N, 0));
    tbb::parallel_for(tbb::blocked_range<int>(0, N), [&](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i != r.end(); ++i) {
            for (int j = 0; j < N; ++j) {
                SD[i][j] = (S[i][j] * diag[j]) % MOD;
            }
        }
        });

    // 5. Final multiplication
    // Note: Ensure your matrix_inverse() and operator* are also using TBB internally!
    return SD * matrix_inverse(S);
}







// Helper to multiply two polynomials over F_p
vector<long long> multiply_polys(const vector<long long>& A, const vector<long long>& B, long long MOD) {
    int degA = A.size() - 1;
    int degB = B.size() - 1;
    vector<long long> C(degA + degB + 1, 0);
    for (int i = 0; i <= degA; ++i) {
        for (int j = 0; j <= degB; ++j) {
            C[i + j] = (C[i + j] + A[i] * B[j]) % MOD;
        }
    }
    return C;
}

// 1. New Struct and Optimized Hessenberg Reduction (Tracks V with Cache-Friendly Transpose)
struct HessenbergResult {
    vector<vector<long long>> H;
    vector<vector<long long>> V;
};

HessenbergResult reduce_to_hessenberg_with_V(vector<vector<long long>> H, long long MOD) {
    int n = H.size();

    // Maintain V as its transpose (V_T) so column operations become fast row operations
    vector<vector<long long>> V_T(n, vector<long long>(n, 0));
    for (int i = 0; i < n; ++i) V_T[i][i] = 1;

    for (int r = 0; r < n - 2; ++r) {
        int pivot = r + 1;
        while (pivot < n && H[pivot][r] == 0) pivot++;
        if (pivot == n) continue;

        if (pivot != r + 1) {
            swap(H[r + 1], H[pivot]);
            for (int i = 0; i < n; ++i)
                swap(H[i][r + 1], H[i][pivot]);
            swap(V_T[r + 1], V_T[pivot]);
        }

        long long inv = inverse(H[r + 1][r]);
        for (int i = r + 2; i < n; ++i) {
            if (H[i][r] != 0) {
                long long factor = (H[i][r] * inv) % MOD;

                // Row operations on H
                for (int j = r; j < n; ++j) {
                    H[i][j] = (H[i][j] - factor * H[r + 1][j]) % MOD;
                    if (H[i][j] < 0) H[i][j] += MOD;
                }

                // Column operations on H
                for (int j = 0; j < n; ++j) {
                    H[j][r + 1] = (H[j][r + 1] + factor * H[j][i]) % MOD;
                }

                // Column operations on V (Using V_T makes this a cache-friendly horizontal scan)
                for (int j = 0; j < n; ++j) {
                    V_T[r + 1][j] = (V_T[r + 1][j] + factor * V_T[i][j]) % MOD;
                }
            }
        }
    }

    // Transpose V_T back to V to return standard row-major format
    vector<vector<long long>> V(n, vector<long long>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            V[i][j] = V_T[j][i];
    return { H, V };
}

// 2. Compute the characteristic polynomial of a single unreduced Hessenberg block using Krylov
vector<long long> get_char_poly_unreduced_hessenberg(const vector<vector<long long>>& H, int start, int end, long long MOD) {
    int m = end - start + 1;
    if (m == 1) {
        // Char poly of 1x1 matrix [h] is x - h
        return { (MOD - H[start][start]) % MOD, 1 };
    }

    // Krylov sequence starting with e_1 = [1, 0, ..., 0]^T
    vector<vector<long long>> K(m, vector<long long>(m, 0));
    vector<long long> curr(m, 0);
    curr[0] = 1;

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) K[i][j] = curr[i];

        // Compute next_v = H_block * curr
        vector<long long> next_v(m, 0);
        for (int i = 0; i < m; ++i) {
            for (int k = 0; k < m; ++k) {
                next_v[i] = (next_v[i] + H[start + i][start + k] * curr[k]) % MOD;
            }
        }
        curr = next_v;
    }

    // Solve K * c = -H^m * e_1
    vector<long long> target(m, 0);
    for (int i = 0; i < m; ++i) target[i] = (MOD - curr[i]) % MOD;

    // Gaussian Elimination (Guaranteed to succeed because block is unreduced Hessenberg!)
    for (int i = 0; i < m; ++i) {
        int pivot = i;
        while (pivot < m && K[pivot][i] == 0) pivot++;

        swap(K[i], K[pivot]);
        swap(target[i], target[pivot]);

        long long inv_pivot = inverse(K[i][i]);
        for (int j = i; j < m; ++j) K[i][j] = (K[i][j] * inv_pivot) % MOD;
        target[i] = (target[i] * inv_pivot) % MOD;

        for (int k = i + 1; k < m; ++k) {
            long long factor = K[k][i];
            for (int j = i; j < m; ++j) {
                K[k][j] = (K[k][j] - factor * K[i][j]) % MOD;
                if (K[k][j] < 0) K[k][j] += MOD;
            }
            target[k] = (target[k] - factor * target[i]) % MOD;
            if (target[k] < 0) target[k] += MOD;
        }
    }

    vector<long long> c(m, 0);
    for (int i = m - 1; i >= 0; --i) {
        c[i] = target[i];
        for (int j = i + 1; j < m; ++j) {
            c[i] = (c[i] - K[i][j] * c[j]) % MOD;
            if (c[i] < 0) c[i] += MOD;
        }
    }

    c.push_back(1); // Explicitly append the leading 1 for x^m
    return c;
}

// 3. Modified Char Poly Wrapper (Accepts pre-reduced H to avoid double-work)
inline vector<long long> get_char_poly_krylov(const vector<vector<long long>>& H, long long MOD) {
    int n = H.size();
    vector<long long> total_poly = { 1 };
    int block_start = 0;

    for (int i = 0; i < n; ++i) {
        if (i == n - 1 || H[i + 1][i] == 0) {
            vector<long long> block_poly = get_char_poly_unreduced_hessenberg(H, block_start, i, MOD);
            total_poly = multiply_polys(total_poly, block_poly, MOD);
            block_start = i + 1;
        }
    }
    return total_poly;
}

// 4. Fast O(n^2) Nullspace for Upper Hessenberg Matrices (With Deferred Modulo)
vector<vector<long long>> get_hessenberg_nullspace(const vector<vector<long long>>& H, long long lambda, long long MOD) {
    int n = H.size();
    if (n == 0) return {};

    vector<vector<long long>> M = H;
    for (int i = 0; i < n; ++i) {
        M[i][i] = (M[i][i] - lambda) % MOD;
        if (M[i][i] < 0) M[i][i] += MOD;
    }

    vector<int> pivot_col_to_row(n, -1);
    int row = 0;

    for (int col = 0; col < n && row < n; ++col) {
        int sel = row;
        while (sel < n && M[sel][col] == 0) sel++;
        if (sel == n) continue;

        if (sel != row) swap(M[row], M[sel]);
        pivot_col_to_row[col] = row;

        long long inv = inverse(M[row][col]);
        for (int j = col; j < n; ++j) M[row][j] = (M[row][j] * inv) % MOD;

        for (int i = row + 1; i < n; ++i) {
            if (M[i][col] != 0) {
                long long factor = M[i][col];
                for (int j = col; j < n; ++j) {
                    long long sub = (factor * M[row][j]) % MOD;
                    M[i][j] = (M[i][j] - sub + MOD) % MOD;
                }
            }
        }
        row++;
    }

    vector<vector<long long>> basis;
    for (int free_col = 0; free_col < n; ++free_col) {
        if (pivot_col_to_row[free_col] == -1) {
            vector<long long> vec(n, 0);
            vec[free_col] = 1;

            // Sequential Back-substitution using Deferred Modulo
            for (int c = free_col - 1; c >= 0; --c) {
                int r = pivot_col_to_row[c];
                if (r != -1) {
                    unsigned long long val = 0;
                    for (int j = c + 1; j <= free_col; ++j) {
                        if (vec[j]) {
                            val += (unsigned long long)M[r][j] * vec[j];
                            if (val >> 61) val %= MOD;
                        }
                    }
                    vec[c] = (MOD - (val % MOD)) % MOD;
                }
            }
            basis.push_back(vec);
        }
    }
    return basis;
}

namespace PolyMath {
    // Helper to strip leading zeros
    void trim(vector<long long>& p) {
        while (p.size() > 1 && p.back() == 0) p.pop_back();
        if (p.empty()) p.push_back(0);
    }

    // O(N^2) Polynomial Multiplication
    vector<long long> mul(const vector<long long>& A, const vector<long long>& B, long long MOD) {
        if (A.empty() || B.empty()) return { 0 };
        vector<long long> C(A.size() + B.size() - 1, 0);
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < B.size(); ++j) {
                C[i + j] = (C[i + j] + A[i] * B[j]) % MOD;
            }
        }
        trim(C);
        return C;
    }

    // O(N^2) Polynomial Long Division (Returns remainder A % B)
    vector<long long> mod(vector<long long> A, const vector<long long>& B, long long MOD) {
        trim(A);
        vector<long long> b = B; trim(b);
        if (b.size() == 1 && b[0] == 0) return { 0 }; // Div by zero fallback
        if (A.size() < b.size()) return A;

        long long inv_lead = inverse(b.back()); // Needs your global inverse()
        for (int i = A.size() - 1; i >= (int)b.size() - 1; --i) {
            if (A[i] == 0) continue;
            long long factor = (A[i] * inv_lead) % MOD;
            for (size_t j = 0; j < b.size(); ++j) {
                long long sub = (factor * b[j]) % MOD;
                A[i - b.size() + 1 + j] = (A[i - b.size() + 1 + j] - sub + MOD) % MOD;
            }
        }
        trim(A);
        return A;
    }

    // O(N^2 * log K) Binary Exponentiation for Polynomials: (base^exp) % P
    vector<long long> pow_mod(vector<long long> base, long long exp, const vector<long long>& P, long long MOD) {
        vector<long long> res = { 1 };
        while (exp > 0) {
            if (exp % 2 == 1) res = mod(mul(res, base, MOD), P, MOD);
            base = mod(mul(base, base, MOD), P, MOD);
            exp /= 2;
        }
        return res;
    }

    // Euclidean Algorithm for Polynomial GCD
    vector<long long> gcd(vector<long long> A, vector<long long> B, long long MOD) {
        trim(A); trim(B);
        while (!(B.size() == 1 && B[0] == 0)) {
            vector<long long> R = mod(A, B, MOD);
            A = B;
            B = R;
        }
        // Normalize so leading coefficient is 1
        if (A.size() > 0 && A.back() != 1) {
            long long inv = inverse(A.back());
            for (long long& x : A) x = (x * inv) % MOD;
        }
        return A;
    }
}

// Internal recursive root finder
void cz_split(const vector<long long>& poly, vector<long long>& roots, long long MOD) {
    if (poly.size() == 2) {
        // Degree 1 polynomial: c_1 * x + c_0 = 0 => x = -c_0 / c_1
        long long root = (MOD - poly[0]) % MOD;
        long long inv_c1 = inverse(poly[1]);
        roots.push_back((root * inv_c1) % MOD);
        return;
    }
    if (poly.size() <= 1) return;

    // Pick a random polynomial a(x) of degree < deg(poly)
    vector<long long> a(poly.size() - 1);
    for (size_t i = 0; i < a.size(); ++i) a[i] = get_rand(MOD);

    // Compute d(x) = gcd(a(x)^((p-1)/2) - 1, poly(x))
    vector<long long> a_pow = PolyMath::pow_mod(a, (MOD - 1) / 2, poly, MOD);
    a_pow[0] = (a_pow[0] - 1 + MOD) % MOD; // Subtract 1

    vector<long long> d = PolyMath::gcd(a_pow, poly, MOD);

    // If it successfully split the polynomial into non-trivial factors, recurse
    if (d.size() > 1 && d.size() < poly.size()) {
        cz_split(d, roots, MOD);
        // poly(x) / d(x) = the other half of the roots
        cz_split(PolyMath::mod(poly, d, MOD), roots, MOD);
    }
    else {
        // Failed to split (50% chance). Try again in the next stack frame.
        cz_split(poly, roots, MOD);
    }
}

// --- The New Optimal Root Finder ---
vector<long long> find_roots_Fp(const vector<long long>& poly, long long MOD) {
    vector<long long> p = poly;
    PolyMath::trim(p);
    if (p.size() <= 1) return {};

    // 1. Compute x^p (mod P) to isolate distinct roots
    vector<long long> X = { 0, 1 }; // The polynomial f(x) = x
    vector<long long> x_pow_p = PolyMath::pow_mod(X, MOD, p, MOD);

    // 2. Subtract x: (x^p - x) mod P
    x_pow_p[0] = (x_pow_p[0] - 0 + MOD) % MOD; // Constant term
    if (x_pow_p.size() < 2) x_pow_p.resize(2, 0);
    x_pow_p[1] = (x_pow_p[1] - 1 + MOD) % MOD; // x term
    PolyMath::trim(x_pow_p);

    // 3. g(x) = gcd(x^p - x, P(x)). This polynomial contains only the distinct linear factors.
    vector<long long> g = PolyMath::gcd(x_pow_p, p, MOD);

    // 4. Split g(x) using Cantor-Zassenhaus
    vector<long long> roots;
    cz_split(g, roots, MOD);
    return roots;
}

// 5. Updated Main Wrapper (Glues everything together in O(n^3) with deferred modulo mapping)
void matrix_diagonalize_krylov(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, long long MOD) {
    int n = A.size();
    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    // Phase 1: Reduce to Hessenberg AND track similarity transformations (O(n^3))
    HessenbergResult hess = reduce_to_hessenberg_with_V(A, MOD);

    // Phase 2: Get characteristic polynomial from reduced H (O(n^3))
    vector<long long> poly = get_char_poly_krylov(hess.H, MOD);

    // Phase 3: Find eigenvalues (Expected O(n^2 log P))
    vector<long long> roots = find_roots_Fp(poly, MOD);

    // Phase 4: Find Eigenvectors of H and map back to A in O(n^3) total
    int col_idx = 0;
    for (long long lambda : roots) {
        // Runs strictly in O(n^2) per eigenvalue
        vector<vector<long long>> H_basis = get_hessenberg_nullspace(hess.H, lambda, MOD);

        for (const auto& y : H_basis) {
            if (col_idx >= n) break;

            // Map back: x = V * y (Takes O(n^2) with deferred modulo)
            for (int i = 0; i < n; ++i) {
                unsigned long long sum = 0;
                for (int j = 0; j < n; ++j) {
                    if (y[j]) {
                        sum += (unsigned long long)hess.V[i][j] * y[j];
                        if (sum >> 61) sum %= MOD;
                    }
                }
                S[i][col_idx] = sum % MOD;
            }
            D[col_idx][col_idx] = lambda;
            col_idx++;
        }
        if (col_idx >= n) break;
    }
}







// --- Experiment B: Full Traditional vs Full Proposed Algorithm ---
inline void EX_B() {
    // 1. Open the text file for writing
    std::ofstream outfile("experiment_B_results.txt");
    if (!outfile.is_open()) {
        printf("Error: Could not open experiment_B_results.txt for writing.\n");
        return;
    }

    // Helper lambda to print to both console and file simultaneously
    auto print_and_write = [&](const std::string& text) {
        printf("%s", text.c_str());
        outfile << text;
        };

    print_and_write("Initializing inverse cache for field F_p...\n");
    for (long long i = 1; i < MOD; ++i)
        inverse(i);
    print_and_write("Cache initialized.\n\n");

    //vector<int> N_values = { 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192 };
    vector<int> N_values = { 512, 1024, 2048, 4096, 8192 };
    //std::vector<int> N_values;
    //for (int i = 1; i < 20; ++i)
    //    N_values.push_back(i * 100);

    print_and_write("=========================================================================================================\n");
    print_and_write("Experiment B: Full Krylov Diagonalization vs Proposed Algorithm (Fixed p)\n");
    print_and_write("=========================================================================================================\n");
    print_and_write("N\tTrials\tTrad Time (s)\tTrad Exp (x)\tProp Time (s)\tProp Exp (x)\tSpeedup\n");
    print_and_write("---------------------------------------------------------------------------------------------------------\n");

    int prev_N = 0;
    double prev_avg_krylov = 0.0;
    double prev_avg_henry = 0.0;

    for (int N : N_values) {
        int num_trials = (N <= 256) ? 10 : (N <= 512) ? 5 : 2;

        double total_time_krylov = 0.0;
        double total_time_henry = 0.0;

        // Ensure both algorithms have fresh output matrices
        std::vector<std::vector<long long>> S_krylov, D_krylov;
        std::vector<std::vector<long long>> S_henry, D_henry;

        for (int trial = 0; trial < num_trials; ++trial) {
            printf("G");
            std::vector<std::vector<long long>> A = generate_test_matrix(N, MOD);

            // --- Time the Full Traditional Baseline ---
            printf("K");
            auto start_krylov = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_krylov(A, S_krylov, D_krylov, MOD);
            auto end_krylov = std::chrono::high_resolution_clock::now();

            if (A * S_krylov != S_krylov * D_krylov) {
                printf("Krylov Diagonalization Error: A*S != S*D\n");
                exit(1);
            }
            bool zero_check = true;
            for (int i = 0; i < N; ++i) {
                if (S_krylov[i][0]) {
                    zero_check = false;
                    break;
                }
            }
            if (zero_check) {
                printf("Krylov Diagonalization Error: S_krylov has zero column\n");
                exit(1);
            }

            std::chrono::duration<double> duration_krylov = end_krylov - start_krylov;
            total_time_krylov += duration_krylov.count();

            // --- Time the Proposed Algorithm ---
            printf("H");
            auto start_henry = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_henry(A, S_henry, D_henry, false);
            auto end_henry = std::chrono::high_resolution_clock::now();

            if (A * S_henry != S_henry * D_henry) {
                printf("Henry Diagonalization Error: A*S != S*D\n");
                exit(1);
            }
            zero_check = true;
            for (int i = 0; i < N; ++i) {
                if (S_henry[i][0]) {
                    zero_check = false;
                    break;
                }
            }
            if (zero_check) {
                printf("Henry Diagonalization Error: S_henry has zero column\n");
                exit(1);
            }

            std::chrono::duration<double> duration_henry = end_henry - start_henry;
            total_time_henry += duration_henry.count();
        }
		printf("\n");

        double avg_krylov = total_time_krylov / num_trials;
        double avg_henry = total_time_henry / num_trials;
        double speedup = avg_krylov / avg_henry;

        // 2. Format the numbers safely into a buffer to retain your exact formatting
        char buffer[256];

        if (prev_N > 0 && prev_avg_krylov > 0.0 && prev_avg_henry > 0.0) {
            // Calculate exponents
            double exp_krylov = std::log(avg_krylov / prev_avg_krylov) / std::log((double)N / prev_N);
            double exp_henry = std::log(avg_henry / prev_avg_henry) / std::log((double)N / prev_N);

            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%.2fx\n",
                N, num_trials, avg_krylov, exp_krylov, avg_henry, exp_henry, speedup);
        }
        else {
            // First run, no previous data
            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\tN/A\t\t%lf\tN/A\t\t%.2fx\n",
                N, num_trials, avg_krylov, avg_henry, speedup);
        }

        // 3. Send formatted buffer to both console and file
        print_and_write(buffer);

        // Store current values for the next iteration
        prev_N = N;
        prev_avg_krylov = avg_krylov;
        prev_avg_henry = avg_henry;
    }

    print_and_write("=========================================================================================================\n");

    // 4. Close the file and notify the user
    outfile.close();
    printf("Results successfully saved to 'experiment_B_results.txt'\n");
}



int main()
{
    //MOD = 1000000007;         //2*500000003         worst distributed
    //MOD = 100000007;          //2*491*101833
    //MOD = 100663291;          //2*3*3*3*5*7*13*16*241
    //MOD = 131071;             //2*3*5*17*257
    //MOD = 524287;             //2*3*3*3*7*19*73     well distributed
    //MOD = 65537;              //2^16
    //MOD = 653659;             //2*3*108943
    //MOD = 257;			    //2*2*2*2*2*2*2*2
    //MOD = 101;                //2*2*5*5

    //Safe primes
    MOD = 565127;


    Initiation();

    EX_B();


    return 0;
}
