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
#include <stdexcept>
#include <exception>

using namespace std;
long long MOD = 10e9; //should be inside of INT32 range. It's GF so MOD must be a prime.
long long primitive;
vector<int> int_inverse;
vector<int> seeds; // primitive ^ seeds[i] = i
vector<int> cubic_z;
vector<long long> MOD_decompose;
vector<long long> MOD_divisors;
vector<vector<int>> ones_roots; // 1^(1/i) = ones_roots.front() ~ back()

inline long long get_rand(long long mod) {
    static std::atomic<int> seed_counter{ 1 };
    thread_local mt19937 rng(
        std::chrono::steady_clock::now().time_since_epoch().count() +
        (seed_counter.fetch_add(1, std::memory_order_relaxed) * 19937)
    );
    uniform_int_distribution<long long> dist(0, mod - 1);
    return dist(rng);
}

std::unique_ptr<tbb::spin_mutex[]> inverse_locks(new tbb::spin_mutex[1024]);
inline long long inverse(long long a) {
    if (!a) {
        throw std::invalid_argument("Integer Inverse Error: 0 has no inverse.");
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
        // Use assign or clear to wipe previous MOD data
        vec.assign(size, 0);
        inFile.read(reinterpret_cast<char*>(vec.data()), size * sizeof(decltype(vec[0])));
        };

    inFile.read(reinterpret_cast<char*>(&primitive), sizeof(primitive));

    readVector(int_inverse);
    readVector(seeds);
    readVector(MOD_decompose);
    readVector(MOD_divisors);
    readVector(cubic_z);

    size_t outerSize;
    inFile.read(reinterpret_cast<char*>(&outerSize), sizeof(outerSize));

    // Clear out the vector of vectors before resizing
    ones_roots.clear();
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
        int_inverse.assign(MOD, 0);
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

        seeds.assign(MOD, 0);
        ones_roots.clear();
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
        cubic_z.assign(MOD, 0);
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

#define TILE_SIZE 16
__global__ void matrixMulModTiledKernel(const long long* A, const long long* B, long long* C, int rowsA, int colsA, int colsB, long long mod) {
    // Allocate Shared Memory for the tile
    __shared__ long long tile_A[TILE_SIZE][TILE_SIZE];
    __shared__ long long tile_B[TILE_SIZE][TILE_SIZE];

    int row = blockIdx.y * TILE_SIZE + threadIdx.y;
    int col = blockIdx.x * TILE_SIZE + threadIdx.x;

    unsigned long long sum = 0;

    // Loop over the tiles of the input matrices
    for (int t = 0; t < (colsA + TILE_SIZE - 1) / TILE_SIZE; ++t) {
        
        // Load data into shared memory (with bounds checking)
        if (row < rowsA && t * TILE_SIZE + threadIdx.x < colsA)
            tile_A[threadIdx.y][threadIdx.x] = A[row * colsA + t * TILE_SIZE + threadIdx.x];
        else
            tile_A[threadIdx.y][threadIdx.x] = 0;

        if (t * TILE_SIZE + threadIdx.y < colsA && col < colsB)
            tile_B[threadIdx.y][threadIdx.x] = B[(t * TILE_SIZE + threadIdx.y) * colsB + col];
        else
            tile_B[threadIdx.y][threadIdx.x] = 0;

        __syncthreads(); // Wait for all threads to finish loading the tile

        // Multiply the tile
        for (int k = 0; k < TILE_SIZE; ++k) {
            sum += (unsigned long long)tile_A[threadIdx.y][k] * tile_B[k][threadIdx.x];
            if (sum >> 61) sum %= mod; // Deferred modulo
        }

        __syncthreads(); // Wait before loading the next tile
    }

    // Write final result to global memory
    if (row < rowsA && col < colsB) {
        C[row * colsB + col] = sum % mod;
    }
}
inline vector<vector<long long>> operator*(const vector<vector<long long>>& a, const vector<vector<long long>>& b) {
    if (a.empty() || b.empty() || a[0].empty() || b[0].empty()) return {};
    int rowsA = a.size(), colsA = a[0].size(), colsB = b[0].size();

    if (colsA != b.size()) {
        throw std::invalid_argument("Matrix Multiplication Error: Matrix size does not match");
    }

    // CPU fallback for small matrices
    if (rowsA <= 20 && colsA <= 20 && colsB <= 20) {
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

    // Flatten input arrays
    vector<long long> flat_A(rowsA * colsA), flat_B(colsA * colsB);
    for (int i = 0; i < rowsA; ++i) copy(a[i].begin(), a[i].end(), flat_A.begin() + i * colsA);
    for (int i = 0; i < colsA; ++i) copy(b[i].begin(), b[i].end(), flat_B.begin() + i * colsB);

    size_t bytesA = rowsA * colsA * sizeof(long long);
    size_t bytesB = colsA * colsB * sizeof(long long);
    size_t bytesC = rowsA * colsB * sizeof(long long);

    // 1. Create stream
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    long long* d_A, * d_B, * d_C;
    // 2. Async memory allocations
    cudaMallocAsync(&d_A, bytesA, stream);
    cudaMallocAsync(&d_B, bytesB, stream);
    cudaMallocAsync(&d_C, bytesC, stream);

    // 3. Async memory copies to Device
    cudaMemcpyAsync(d_A, flat_A.data(), bytesA, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_B, flat_B.data(), bytesB, cudaMemcpyHostToDevice, stream);

    dim3 threads(16, 16);
    dim3 blocks((colsB + 15) / 16, (rowsA + 15) / 16);

    // 4. Launch kernel into stream
    matrixMulModTiledKernel << <blocks, threads, 0, stream >> > (d_A, d_B, d_C, rowsA, colsA, colsB, MOD);

    // 5. Async memory copy back to Host, then strictly sync ONLY this stream
    vector<long long> flat_C(rowsA * colsB);
    cudaMemcpyAsync(flat_C.data(), d_C, bytesC, cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    // Async memory free and stream destruction
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_B, stream);
    cudaFreeAsync(d_C, stream);
    cudaStreamDestroy(stream);

    // Reconstruct 2D Vector
    vector<vector<long long>> R(rowsA, vector<long long>(colsB));
    for (int i = 0; i < rowsA; ++i) {
        copy(flat_C.begin() + i * colsB, flat_C.begin() + (i + 1) * colsB, R[i].begin());
    }

    return R;
}

#define TILE_SIZE2 256
__global__ void matrixVectorMulModKernel(const long long* __restrict__ A, const long long* __restrict__ B, long long* __restrict__ C, int rowsA, int colsA, long long mod) {
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long sum = 0;

    // Loop over vector B in TILE_SIZE chunks
    for (int t = 0; t < (colsA + TILE_SIZE2 - 1) / TILE_SIZE2; ++t) {

        // Allocate shared memory for the chunk of vector B
        __shared__ long long shared_B[TILE_SIZE2];

        // Let the threads cooperatively load a chunk of B into shared memory
        int b_idx = t * TILE_SIZE2 + threadIdx.x;
        if (b_idx < colsA) {
            shared_B[threadIdx.x] = B[b_idx];
        }

        // Wait for all threads in the block to finish loading shared_B
        __syncthreads();

        // Compute the partial dot product using the ultra-fast shared memory
        if (row < rowsA) {
            int limit = (t + 1) * TILE_SIZE2 > colsA ? colsA - t * TILE_SIZE2 : TILE_SIZE2;
            for (int j = 0; j < limit; ++j) {
                sum += (unsigned long long)A[row * colsA + t * TILE_SIZE2 + j] * shared_B[j];
                if (sum >> 61) sum %= mod; // Deferred modulo trick
            }
        }

        // Wait for all threads to finish computing before overwriting shared_B in the next loop
        __syncthreads();
    }

    // Write final result to global memory
    if (row < rowsA) {
        C[row] = sum % mod;
    }
}
inline vector<long long> operator*(const vector<vector<long long>>& a, const vector<long long>& b) {
    if (a.empty() || a[0].empty() || b.empty()) return {};
    if (a[0].size() != b.size()) {
        throw std::invalid_argument("Matrix Vector Multiplication Error: Matrix and Vector's size do not match");
    }

    int rowsA = a.size(), colsA = a[0].size();

    // CPU fallback for small workloads
    if (rowsA <= 20 && colsA <= 20) {
        vector<long long> R(rowsA, 0);
        for (int i = 0; i < rowsA; ++i) {
            unsigned long long sum = 0;
            for (int j = 0; j < colsA; ++j) {
                sum += (unsigned long long)a[i][j] * b[j];
                // Fixed: Apply the fast deferred modulo trick to the CPU as well!
                if (sum >> 61) sum %= MOD;
            }
            R[i] = sum % MOD;
        }
        return R;
    }

    // Flatten A
    vector<long long> flat_A(rowsA * colsA), R(rowsA);
    for (int i = 0; i < rowsA; ++i) {
        copy(a[i].begin(), a[i].end(), flat_A.begin() + i * colsA);
    }

    size_t bytesA = rowsA * colsA * sizeof(long long);
    size_t bytesB = colsA * sizeof(long long);
    size_t bytesC = rowsA * sizeof(long long);

    long long* d_A, * d_B, * d_C;

    // Dedicated stream for TBB concurrency
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    // Async memory allocations
    cudaMallocAsync(&d_A, bytesA, stream);
    cudaMallocAsync(&d_B, bytesB, stream);
    cudaMallocAsync(&d_C, bytesC, stream);

    // Async transfers to device
    cudaMemcpyAsync(d_A, flat_A.data(), bytesA, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_B, b.data(), bytesB, cudaMemcpyHostToDevice, stream);

    // Launch tiled kernel into the stream
    int threads = TILE_SIZE; // Matches TILE_SIZE for optimal shared memory loading
    int blocks = (rowsA + threads - 1) / threads;
    matrixVectorMulModKernel << <blocks, threads, 0, stream >> > (d_A, d_B, d_C, rowsA, colsA, MOD);

    // Async transfer back to host
    cudaMemcpyAsync(R.data(), d_C, bytesC, cudaMemcpyDeviceToHost, stream);

    // CPU waits ONLY for this specific stream to finish
    cudaStreamSynchronize(stream);

    // Async cleanup
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_B, stream);
    cudaFreeAsync(d_C, stream);
    cudaStreamDestroy(stream);

    return R;
}
inline long long operator * (const vector<long long>& a, const vector<long long>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector Dot Product Error: Vector size does not match");
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
        throw std::invalid_argument("Matrix Addition Error: Matrix size does not match");
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
        throw std::invalid_argument("Matrix Subtraction Error: Matrix size does not match");
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
        throw std::invalid_argument("Vector Addition Error: Vector size does not match");
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
        throw std::invalid_argument("Vector Subtraction Error: Vector size does not match");
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
__global__ void initIdentityKernel(long long* I, int n) {
    int r = blockIdx.y * blockDim.y + threadIdx.y;
    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (r < n && c < n) I[r * n + c] = (r == c) ? 1 : 0;
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
inline vector<vector<long long>> matrix_power_gpu(const vector<vector<long long>>& a, unsigned long long p) {
    if (a.empty() || a[0].empty()) return {};
    int N = a.size();

    // CPU fallback for small matrices (avoids GPU launch latency)
    if (N <= 20) {
        vector<vector<long long>> res = I_n(N); // Your existing identity function
        vector<vector<long long>> base = a;
        while (p) {
            if (p & 1) res = res * base; // Uses the CPU fallback inside operator*
            p >>= 1;
            if (!p) break;
            base = base * base;
        }
        return res;
    }

    // Flatten input matrix
    vector<long long> flat_a(N * N);
    for (int i = 0; i < N; ++i) {
        copy(a[i].begin(), a[i].end(), flat_a.begin() + i * N);
    }

    size_t bytes = N * N * sizeof(long long);

    // Dedicated stream for TBB concurrency
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    long long* d_A, * d_res, * d_temp;
    cudaMallocAsync(&d_A, bytes, stream);
    cudaMallocAsync(&d_res, bytes, stream);
    cudaMallocAsync(&d_temp, bytes, stream); // Workspace buffer

    // Copy base matrix to device
    cudaMemcpyAsync(d_A, flat_a.data(), bytes, cudaMemcpyHostToDevice, stream);

    // Initialize d_res as Identity Matrix directly on the GPU
    dim3 threads(16, 16);
    dim3 blocks((N + 15) / 16, (N + 15) / 16);
    initIdentityKernel << <blocks, threads, 0, stream >> > (d_res, N);

    // ---------------------------------------------------------
    // THE KERNEL LOOP (Zero PCIe Transfers here!)
    // ---------------------------------------------------------
    while (p) {
        if (p & 1) {
            // d_temp = d_res * d_A
            matrixMulModTiledKernel << <blocks, threads, 0, stream >> > (d_res, d_A, d_temp, N, N, N, MOD);
            // Swap pointers: d_res now holds the new result, d_temp becomes the free buffer
            std::swap(d_res, d_temp);
        }
        p >>= 1;
        if (!p) break;

        // d_temp = d_A * d_A
        matrixMulModTiledKernel << <blocks, threads, 0, stream >> > (d_A, d_A, d_temp, N, N, N, MOD);
        // Swap pointers: d_A now holds the squared matrix, d_temp becomes the free buffer
        std::swap(d_A, d_temp);
    }

    // ---------------------------------------------------------
    // Retrieve Final Result
    // ---------------------------------------------------------
    vector<long long> flat_res(N * N);

    // We only copy memory BACK to the CPU when the entire exponentiation is done
    cudaMemcpyAsync(flat_res.data(), d_res, bytes, cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream); // CPU waits here

    // Async cleanup
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_res, stream);
    cudaFreeAsync(d_temp, stream);
    cudaStreamDestroy(stream);

    // Reconstruct 2D matrix
    vector<vector<long long>> res(N, vector<long long>(N));
    for (int i = 0; i < N; ++i) {
        copy(flat_res.begin() + i * N, flat_res.begin() + (i + 1) * N, res[i].begin());
    }

    return res;
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
            throw std::invalid_argument("Matrix Chop Error: Block bounds exceed matrix dimensions");
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
            throw std::invalid_argument("Matrix Partial Multiply Error: Block bounds exceed matrix dimensions");
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





// Device-side modular inverse
__device__ long long inverse_gpu(long long a, long long m) {
    long long m0 = m, t, q;
    long long x0 = 0, x1 = 1;
    if (m == 1) return 0;
    while (a > 1) {
        q = a / m; t = m; m = a % m; a = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    return x1 < 0 ? x1 + m0 : x1;
}

// ------------------------------------------------------------------
// 1. Pivot Search 
// ------------------------------------------------------------------
__global__ void findPivotKernel(const long long* A, int m, int cols_A, int static_row, const int* d_dynamic_row, int col, int* out_pivot_row, long long* out_pivot_val, const int* d_info) {
    if (d_info != nullptr && *d_info != 0) return;

    // Use dynamic row if provided (Rank/NullSpace), otherwise static (Inverse/Det)
    int current_row = (d_dynamic_row != nullptr) ? *d_dynamic_row : static_row;

    __shared__ int min_row;
    if (threadIdx.x == 0) min_row = m;
    __syncthreads();

    if (current_row < m) {
        for (int r = current_row + threadIdx.x; r < m; r += blockDim.x) {
            if (A[r * cols_A + col] != 0) {
                atomicMin(&min_row, r);
            }
        }
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        *out_pivot_row = min_row;
        if (out_pivot_val != nullptr) {
            if (min_row < m) *out_pivot_val = A[min_row * cols_A + col];
            else *out_pivot_val = 0;
        }
    }
}

// ------------------------------------------------------------------
// 2. Row Swap 
// ------------------------------------------------------------------
template <bool HAS_AUG>
__global__ void swapRowsKernel(long long* A, long long* B, int m, int cols_A, int cols_B, int static_row, const int* d_dynamic_row, const int* d_pivot_row, const int* d_info) {
    if (d_info != nullptr && *d_info != 0) return;

    int pivot_row = *d_pivot_row;
    if (pivot_row == m) return; // No pivot found in this column, skip

    int current_row = (d_dynamic_row != nullptr) ? *d_dynamic_row : static_row;
    if (current_row == pivot_row) return; // Already in place

    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c < cols_A) {
        long long t = A[current_row * cols_A + c];
        A[current_row * cols_A + c] = A[pivot_row * cols_A + c];
        A[pivot_row * cols_A + c] = t;
    }
    if (HAS_AUG && c < cols_B) {
        long long t = B[current_row * cols_B + c];
        B[current_row * cols_B + c] = B[pivot_row * cols_B + c];
        B[pivot_row * cols_B + c] = t;
    }
}

// ------------------------------------------------------------------
// 3. Normalize Pivot Row 
// ------------------------------------------------------------------
template <bool HAS_AUG>
__global__ void normalizePivotRowKernel(long long* A, long long* B, int m, int cols_A, int cols_B, int static_row, const int* d_dynamic_row, int pivot_col, const int* d_pivot_row, long long mod, const long long* d_pivot_val, const int* d_info) {
    if (d_info != nullptr && *d_info != 0) return;
    if (*d_pivot_row == m) return;

    int current_row = (d_dynamic_row != nullptr) ? *d_dynamic_row : static_row;

    __shared__ long long inv;
    if (threadIdx.x == 0) {
        // [FIX]: Read the safely isolated pivot value. No read-after-write hazard!
        inv = inverse_gpu(*d_pivot_val, mod);
    }
    __syncthreads();

    int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c < cols_A) {
        A[current_row * cols_A + c] = (A[current_row * cols_A + c] * inv) % mod;
    }
    if (HAS_AUG && c < cols_B) {
        B[current_row * cols_B + c] = (B[current_row * cols_B + c] * inv) % mod;
    }
}

// ------------------------------------------------------------------
// 4. Compute Multipliers 
// ------------------------------------------------------------------
__global__ void computeMultipliersKernel(const long long* A, int m, int cols_A, int static_row, const int* d_dynamic_row, int pivot_col, const int* d_pivot_row, long long mod, long long* multipliers, bool both_directions, const int* d_info) {
    if (d_info != nullptr && *d_info != 0) return;
    if (*d_pivot_row == m) return; // No pivot found, skip

    int current_row = (d_dynamic_row != nullptr) ? *d_dynamic_row : static_row;

    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < m) {
        if (j == current_row || (!both_directions && j < current_row)) {
            multipliers[j] = 0;
        }
        else {
            long long val = A[j * cols_A + pivot_col];
            multipliers[j] = (val == 0) ? 0 : (mod - val) % mod;
        }
    }
}

// ------------------------------------------------------------------
// 5. Eliminate Rows 
// ------------------------------------------------------------------
template <bool HAS_AUG>
__global__ void eliminateRowsKernel(long long* A, long long* B, int m, int cols_A, int cols_B, int static_row, const int* d_dynamic_row, const int* d_pivot_row, const long long* multipliers, long long mod, const int* d_info) {
    if (d_info != nullptr && *d_info != 0) return;
    if (*d_pivot_row == m) return; // No pivot found, skip

    int current_row = (d_dynamic_row != nullptr) ? *d_dynamic_row : static_row;

    int r = blockIdx.y * blockDim.y + threadIdx.y;
    int c = blockIdx.x * blockDim.x + threadIdx.x;

    if (r < m) {
        long long mul = multipliers[r];
        if (mul != 0) {
            if (c < cols_A) {
                A[r * cols_A + c] = (A[r * cols_A + c] + A[current_row * cols_A + c] * mul) % mod;
            }
            if (HAS_AUG && c < cols_B) {
                B[r * cols_B + c] = (B[r * cols_B + c] + B[current_row * cols_B + c] * mul) % mod;
            }
        }
    }
}

// ------------------------------------------------------------------
// HELPER: Update State for Determinant / Inverse
// ------------------------------------------------------------------
__global__ void updateDetInfoKernel(const int* d_pivot_row, const long long* d_pivot_val, int static_row, int m, long long* d_det, int* d_info, long long mod) {
    if (d_info != nullptr && *d_info != 0) return;
    if (threadIdx.x != 0 || blockIdx.x != 0) return;

    if (*d_pivot_row == m) {
        if (d_info != nullptr) *d_info = 1; // Flag Singular
        if (d_det != nullptr) *d_det = 0;   // Det is 0
    }
    else if (d_det != nullptr) {
        if (*d_pivot_row != static_row) {
            *d_det = (mod - *d_det) % mod;  // Swap flips sign
        }
        *d_det = (*d_det * (*d_pivot_val)) % mod;
    }
}

// ------------------------------------------------------------------
// HELPER: Update State for Rank / Null Space
// ------------------------------------------------------------------
__global__ void updateRankStateKernel(const int* d_pivot_row, int m, int col, int* d_dynamic_row, int* d_piv_array) {
    // Only increment row and record pivot if a valid pivot was found
    if (threadIdx.x != 0 || blockIdx.x != 0) return;
    if (*d_pivot_row < m) {
        if (d_piv_array != nullptr) {
            d_piv_array[*d_dynamic_row] = col;
        }
        (*d_dynamic_row)++;
    }
}




inline size_t matrix_rank_gpu(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D[0].empty()) return 0;
    int m = A_2D.size(), n = A_2D[0].size();

    vector<long long> flat_A(m * n);
    for (int i = 0; i < m; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), flat_A.begin() + i * n);
    }

    cudaStream_t stream;
    cudaStreamCreate(&stream); // Independent stream

    long long* d_A, * d_mults, * d_pivot_val; // [FIX]: Added d_pivot_val
    int* d_pivot_row, * d_row;
    int* h_rank;

    // Async memory allocations
    cudaMallocAsync(&d_A, m * n * sizeof(long long), stream);
    cudaMallocAsync(&d_mults, m * sizeof(long long), stream);
    cudaMallocAsync(&d_pivot_row, sizeof(int), stream);
    cudaMallocAsync(&d_pivot_val, sizeof(long long), stream); // [FIX]: Allocate pivot value isolator

    // NEW: Device-side row counter to break CPU dependency
    cudaMallocAsync(&d_row, sizeof(int), stream);
    cudaMallocHost(&h_rank, sizeof(int)); // Pinned memory

    // Initialize state
    cudaMemsetAsync(d_row, 0, sizeof(int), stream);
    cudaMemcpyAsync(d_A, flat_A.data(), m * n * sizeof(long long), cudaMemcpyHostToDevice, stream);

    int threads1D = 256;
    int blocks1D = (n + threads1D - 1) / threads1D;
    int blocksMult = (m + threads1D - 1) / threads1D;
    dim3 threads2D(16, 16), blocks2D((n + 15) / 16, (m + 15) / 16);

    // ENTIRE LOOP QUEUED ASYNCHRONOUSLY
    // CPU runs through this in microseconds. If d_row reaches m, 
    // the kernels gracefully short-circuit without crashing.
    for (int col = 0; col < n; ++col) {
        // [FIX]: Reset pivot row to -1 to prevent phantom pivots on free columns
        cudaMemsetAsync(d_pivot_row, -1, sizeof(int), stream);

        // 1. Find pivot using dynamic d_row, saving value to d_pivot_val
        findPivotKernel << <1, 256, 0, stream >> > (
            d_A, m, n, 0, d_row, col, d_pivot_row, d_pivot_val, nullptr);

        // 2. Math kernels conditionally execute based on d_pivot_row
        swapRowsKernel<false> << <blocks1D, threads1D, 0, stream >> > (
            d_A, nullptr, m, n, 0, 0, d_row, d_pivot_row, nullptr);

        // [FIX]: Pass d_pivot_val so all blocks read the original uncorrupted value
        normalizePivotRowKernel<false> << <blocks1D, threads1D, 0, stream >> > (
            d_A, nullptr, m, n, 0, 0, d_row, col, d_pivot_row, MOD, d_pivot_val, nullptr);

        computeMultipliersKernel << <blocksMult, threads1D, 0, stream >> > (
            d_A, m, n, 0, d_row, col, d_pivot_row, MOD, d_mults, false, nullptr);

        eliminateRowsKernel<false> << <blocks2D, threads2D, 0, stream >> > (
            d_A, nullptr, m, n, 0, 0, d_row, d_pivot_row, d_mults, MOD, nullptr);

        // 3. Increment d_row AT THE END of the iteration if a pivot was successfully found
        updateRankStateKernel << <1, 1, 0, stream >> > (
            d_pivot_row, m, col, d_row, nullptr);
    }

    // ONLY SYNC ONCE AT THE VERY END to retrieve the rank (which is exactly the final d_row value)
    cudaMemcpyAsync(h_rank, d_row, sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    size_t rank = *h_rank;

    // Cleanup asynchronously
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_mults, stream);
    cudaFreeAsync(d_pivot_row, stream);
    cudaFreeAsync(d_pivot_val, stream); // [FIX]: Free isolated pivot value
    cudaFreeAsync(d_row, stream);
    cudaFreeHost(h_rank);
    cudaStreamDestroy(stream);

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
inline vector<vector<long long>> matrix_inverse_gpu(const vector<vector<long long>>& A_2D) {
    int n = A_2D.size();
    size_t bytes = n * n * sizeof(long long);

    // Flat arrays
    vector<long long> flat_A(n * n), flat_I(n * n);
    for (int i = 0; i < n; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), flat_A.begin() + i * n);
    }

    // 1. Create independent stream for this TBB thread
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    long long* d_A, * d_I, * d_mults, * d_pivot_val; // [FIX]: Added d_pivot_val
    int* d_pivot, * d_info;

    // 2. Async Allocations (No thread blocking, extremely fast)
    cudaMallocAsync(&d_A, bytes, stream);
    cudaMallocAsync(&d_I, bytes, stream);
    cudaMallocAsync(&d_mults, n * sizeof(long long), stream);
    cudaMallocAsync(&d_pivot, sizeof(int), stream);
    cudaMallocAsync(&d_info, sizeof(int), stream);
    cudaMallocAsync(&d_pivot_val, sizeof(long long), stream); // [FIX]: Allocate d_pivot_val

    cudaMemcpyAsync(d_A, flat_A.data(), bytes, cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(d_info, 0, sizeof(int), stream);

    // Concise thread configs
    dim3 b1D((n + 255) / 256), t1D(256);
    dim3 b2D((n + 15) / 16, (n + 15) / 16), t2D(16, 16);

    initIdentityKernel << <b2D, t2D, 0, stream >> > (d_I, n);

    // 3. ENTIRE PIPELINE QUEUED ASYNCHRONOUSLY
    for (int p = 0; p < n; ++p) {
        // [FIX]: Pass d_pivot_val to capture the pivot mathematically safely
        findPivotKernel << <1, 256, 0, stream >> > (d_A, n, n, p, nullptr, p, d_pivot, d_pivot_val, d_info);

        // Optional: Pass it to DetInfoKernel as well for consistency
        updateDetInfoKernel << <1, 1, 0, stream >> > (d_pivot, d_pivot_val, p, n, nullptr, d_info, MOD);

        swapRowsKernel<true> << <b1D, t1D, 0, stream >> > (d_A, d_I, n, n, n, p, nullptr, d_pivot, d_info);

        // [FIX]: Pass d_pivot_val into the normalizer to prevent read-after-write block hazards
        normalizePivotRowKernel<true> << <b1D, t1D, 0, stream >> > (d_A, d_I, n, n, n, p, nullptr, p, d_pivot, MOD, d_pivot_val, d_info);

        computeMultipliersKernel << <b1D, t1D, 0, stream >> > (d_A, n, n, p, nullptr, p, d_pivot, MOD, d_mults, true, d_info);
        eliminateRowsKernel<true> << <b2D, t2D, 0, stream >> > (d_A, d_I, n, n, n, p, nullptr, d_pivot, d_mults, MOD, d_info);
    }

    // Standard stack variable, no expensive OS allocations needed
    int h_info_val = 0;

    // 4. Fetch the result and the error flag
    // (Copying to pageable memory blocks the host, which is fine since we sync immediately after)
    cudaMemcpyAsync(&h_info_val, d_info, sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(flat_I.data(), d_I, bytes, cudaMemcpyDeviceToHost, stream);

    // 5. The ONLY synchronization point
    cudaStreamSynchronize(stream);

    // 6. Queue Async Cleanups BEFORE throwing any exceptions
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_I, stream);
    cudaFreeAsync(d_mults, stream);
    cudaFreeAsync(d_pivot, stream);
    cudaFreeAsync(d_info, stream);
    cudaFreeAsync(d_pivot_val, stream); // [FIX]: Cleanup new pointer
    cudaStreamDestroy(stream);

    // 7. Safe Error Check
    if (h_info_val == 1) {
        printf("Matrix is singular\n");
        throw std::invalid_argument("Matrix is singular");
    }

    // Construct final output
    vector<vector<long long>> I_out(n, vector<long long>(n));
    for (int i = 0; i < n; ++i) {
        copy(flat_I.begin() + i * n, flat_I.begin() + (i + 1) * n, I_out[i].begin());
    }

    return I_out;
}
inline vector<vector<long long>> matrix_inverse(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D.size() != A_2D.front().size()) {
        throw std::invalid_argument("Matrix Inversion Error: Matrix is not square");
    }

    size_t n = A_2D.size();
    vector<long long> A(n * n);
    vector<long long> I(n * n, 0);

    // 1. PARALLELIZE INITIALIZATION
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                copy(A_2D[i].begin(), A_2D[i].end(), A.begin() + i * n);
                I[i * n + i] = 1;
            }
        });

    // 2. FUSE FORWARD AND BACKWARD ELIMINATION
    for (size_t p = 0; p < n; ++p) {
        // Pivot search (Sequential, as we break early on success)
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
                throw std::invalid_argument("Matrix Inversion Error: Matrix is singular");
            }
        }

        long long inv_pivot = inverse(A[p * n + p]);

        // Eliminate ALL other rows (both above and below the pivot) simultaneously
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t j = r.begin(); j != r.end(); ++j) {
                    if (j == p) continue; // Skip the pivot row itself

                    long long val = A[j * n + p];
                    if (val == 0) continue; // Skip if already zero (saves math operations)

                    long long mul = (MOD - val) * inv_pivot % MOD;

                    for (size_t k = p; k < n; ++k)
                        A[j * n + k] = (A[j * n + k] + A[p * n + k] * mul) % MOD;

                    for (size_t k = 0; k < n; ++k)
                        I[j * n + k] = (I[j * n + k] + I[p * n + k] * mul) % MOD;
                }
            });
    }

    // 3. DIAGONAL NORMALIZATION (Extract output directly)
    vector<vector<long long>> I_out(n, vector<long long>(n));
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                long long t = inverse(A[i * n + i]); // The diagonal element
                for (size_t j = 0; j < n; ++j) {
                    I_out[i][j] = (I[i * n + j] * t) % MOD;
                }
            }
        });

    return I_out;
}
inline long long matrix_determinant_gpu(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D.size() != A_2D.front().size()) {
        throw std::invalid_argument("Matrix determinant Error: Matrix is not square");
    }

    int n = A_2D.size();
    vector<long long> flat_A(n * n);
    for (int i = 0; i < n; ++i) {
        copy(A_2D[i].begin(), A_2D[i].end(), flat_A.begin() + i * n);
    }

    // 1. Create a dedicated Stream for this thread
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    long long* d_A, * d_mults, * d_pivot_val, * d_det;
    int* d_pivot_row, * d_info;

    // 2. Async Device Allocations (Non-blocking)
    cudaMallocAsync(&d_A, n * n * sizeof(long long), stream);
    cudaMallocAsync(&d_mults, n * sizeof(long long), stream);
    cudaMallocAsync(&d_pivot_row, sizeof(int), stream);
    cudaMallocAsync(&d_pivot_val, sizeof(long long), stream);

    // Allocate state trackers on the device
    cudaMallocAsync(&d_det, sizeof(long long), stream);
    cudaMallocAsync(&d_info, sizeof(int), stream);

    // Initialize state trackers
    long long h_det_init = 1;
    cudaMemcpyAsync(d_det, &h_det_init, sizeof(long long), cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(d_info, 0, sizeof(int), stream);

    // 3. Async Copy to Device
    cudaMemcpyAsync(d_A, flat_A.data(), n * n * sizeof(long long), cudaMemcpyHostToDevice, stream);

    int t1D = 256, b1D = (n + t1D - 1) / t1D;
    dim3 t2D(16, 16), b2D((n + 15) / 16, (n + 15) / 16);

    // 4. QUEUE ALL KERNELS ASYNCHRONOUSLY
    for (int p = 0; p < n; ++p) {
        // Find the pivot and write it to device pointers
        findPivotKernel << <1, 256, 0, stream >> > (d_A, n, n, p, nullptr, p, d_pivot_row, d_pivot_val, d_info);

        // Update determinant and singularity flag instantly on the GPU
        updateDetInfoKernel << <1, 1, 0, stream >> > (d_pivot_row, d_pivot_val, p, n, d_det, d_info, MOD);

        // Swap rows if needed (Short-circuits if d_info == 1 or d_pivot_row == m)
        swapRowsKernel<false> << <b1D, t1D, 0, stream >> > (d_A, nullptr, n, n, 0, p, nullptr, d_pivot_row, d_info);

        // CPU can still safely evaluate loop boundaries since 'p' and 'n' are strictly CPU constants
        if (p == n - 1) break;

        // Normalize and eliminate
        // [FIX]: Added d_pivot_val before d_info to match the updated kernel signature!
        normalizePivotRowKernel<false> << <b1D, t1D, 0, stream >> > (d_A, nullptr, n, n, 0, p, nullptr, p, d_pivot_row, MOD, d_pivot_val, d_info);

        computeMultipliersKernel << <b1D, t1D, 0, stream >> > (d_A, n, n, p, nullptr, p, d_pivot_row, MOD, d_mults, false, d_info);
        eliminateRowsKernel<false> << <b2D, t2D, 0, stream >> > (d_A, nullptr, n, n, 0, p, nullptr, d_pivot_row, d_mults, MOD, d_info);
    }

    // 5. One single synchronization point at the end
    long long final_det = 0;
    cudaMemcpyAsync(&final_det, d_det, sizeof(long long), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    // 6. Async Device Cleanup
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_mults, stream);
    cudaFreeAsync(d_pivot_row, stream);
    cudaFreeAsync(d_pivot_val, stream);
    cudaFreeAsync(d_det, stream);
    cudaFreeAsync(d_info, stream);

    cudaStreamDestroy(stream);

    return final_det;
}
inline long long matrix_determinant(const vector<vector<long long>>& A_2D) {
    if (A_2D.empty() || A_2D.size() != A_2D.front().size()) {
        throw std::invalid_argument("Matrix determinant Error: Matrix is not square");
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
inline vector<vector<long long>> Null_Space_gpu(const vector<vector<long long>>& A_2D, bool Orth) {
    if (A_2D.empty() || A_2D[0].empty()) return {};

    int m = A_2D.size(), n = A_2D[0].size();
    vector<long long> flat_A(m * n);
    for (int i = 0; i < m; ++i) copy(A_2D[i].begin(), A_2D[i].end(), flat_A.begin() + i * n);

    // ---------------------------------------------------------
    // PHASE 1: GPU RREF (100% Asynchronous Pipeline)
    // ---------------------------------------------------------
    cudaStream_t stream;
    cudaStreamCreate(&stream);

    long long* d_A, * d_mults, * d_pivot_val;
    int* d_pivot_row, * d_dynamic_row, * d_piv_array;

    // Use standard stack/heap allocations instead of slow cudaMallocHost
    int h_rank = 0;
    vector<int> h_piv_array(n, 0);

    // Async Device Allocations
    cudaMallocAsync(&d_A, m * n * sizeof(long long), stream);
    cudaMallocAsync(&d_mults, m * sizeof(long long), stream);
    cudaMallocAsync(&d_pivot_row, sizeof(int), stream);
    cudaMallocAsync(&d_dynamic_row, sizeof(int), stream);
    cudaMallocAsync(&d_piv_array, n * sizeof(int), stream);
    cudaMallocAsync(&d_pivot_val, sizeof(long long), stream); // [FIX]: Isolate pivot value

    // Initial copies and setups
    cudaMemcpyAsync(d_A, flat_A.data(), m * n * sizeof(long long), cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(d_dynamic_row, 0, sizeof(int), stream);

    int threads = 256;
    dim3 threads2D(16, 16), blocks2D((n + 15) / 16, (m + 15) / 16);
    dim3 blocks1D_cols((n + threads - 1) / threads);
    dim3 blocks1D_rows((m + threads - 1) / threads);

    // Enqueue kernels asynchronously
    for (int col = 0; col < n; ++col) {
        // [FIX]: Prevent "Phantom Pivots" by resetting to -1 before searching
        cudaMemsetAsync(d_pivot_row, -1, sizeof(int), stream);

        // [FIX]: Pass d_pivot_val to store the original pivot value safely
        findPivotKernel << <1, 256, 0, stream >> > (d_A, m, n, 0, d_dynamic_row, col, d_pivot_row, d_pivot_val, nullptr);

        swapRowsKernel<false> << <blocks1D_cols, threads, 0, stream >> > (d_A, nullptr, m, n, 0, 0, d_dynamic_row, d_pivot_row, nullptr);

        // [FIX]: Pass d_pivot_val so all blocks read the original uncorrupted value
        normalizePivotRowKernel<false> << <blocks1D_cols, threads, 0, stream >> > (d_A, nullptr, m, n, 0, 0, d_dynamic_row, col, d_pivot_row, MOD, d_pivot_val, nullptr);

        computeMultipliersKernel << <blocks1D_rows, threads, 0, stream >> > (d_A, m, n, 0, d_dynamic_row, col, d_pivot_row, MOD, d_mults, true, nullptr);

        eliminateRowsKernel<false> << <blocks2D, threads2D, 0, stream >> > (d_A, nullptr, m, n, 0, 0, d_dynamic_row, d_pivot_row, d_mults, MOD, nullptr);

        updateRankStateKernel << <1, 1, 0, stream >> > (d_pivot_row, m, col, d_dynamic_row, d_piv_array);
    }

    // Copy back to standard host memory and synchronize
    cudaMemcpyAsync(flat_A.data(), d_A, m * n * sizeof(long long), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(&h_rank, d_dynamic_row, sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(h_piv_array.data(), d_piv_array, n * sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    // Cleanup Async
    cudaFreeAsync(d_A, stream);
    cudaFreeAsync(d_mults, stream);
    cudaFreeAsync(d_pivot_row, stream);
    cudaFreeAsync(d_dynamic_row, stream);
    cudaFreeAsync(d_piv_array, stream);
    cudaFreeAsync(d_pivot_val, stream); // [FIX]: Clean up new allocation
    cudaStreamDestroy(stream);

    // ---------------------------------------------------------
    // PHASE 2: CPU Null Space Extraction & Orthogonalization
    // ---------------------------------------------------------
    int rank = h_rank;
    if (rank == n) return {};

    int null_dim = n - rank;
    vector<vector<long long>> NS(null_dim, vector<long long>(n, 0));
    vector<bool> is_pivot(n, false);
    for (int i = 0; i < rank; i++) is_pivot[h_piv_array[i]] = true;

    for (int col = 0, free_idx = 0; col < n; ++col) {
        if (!is_pivot[col]) {
            NS[free_idx][col] = 1;
            for (int i = 0; i < rank; ++i)
                NS[free_idx][h_piv_array[i]] = (MOD - flat_A[i * n + col]) % MOD;
            free_idx++;
        }
    }

    if (Orth) {
        vector<long long> DP(null_dim, 0);
        for (int i = 0; i < null_dim; ++i) {
            for (int j = 0; j < i; ++j) {
                // OPTIMIZATION: Conditional modulo avoids __int128 and defers expensive % ops
                unsigned long long dot_accum = 0;
                for (int k = 0; k < n; ++k) {
                    dot_accum += (unsigned long long)NS[i][k] * NS[j][k];
                    if (dot_accum >> 61) dot_accum %= MOD;
                }
                long long dot = (long long)(dot_accum % MOD);

                long long c = (dot * inverse(DP[j])) % MOD;
                for (int k = 0; k < n; ++k) {
                    NS[i][k] = (NS[i][k] - ((c * NS[j][k]) % MOD) + MOD) % MOD;
                }
            }

            unsigned long long norm_sq_accum = 0;
            for (int k = 0; k < n; ++k) {
                norm_sq_accum += (unsigned long long)NS[i][k] * NS[i][k];
                if (norm_sq_accum >> 61) norm_sq_accum %= MOD;
            }

            if ((DP[i] = (long long)(norm_sq_accum % MOD)) == 0) {
                throw std::invalid_argument("NullSpace's G-S Process Error: Isotropic vector encountered");
            }
        }
    }

    // OPTIMIZATION: Swap loop dimensions for better cache locality (Write linearly into NS_col[j])
    vector<vector<long long>> NS_col(n, vector<long long>(null_dim));
    tbb::parallel_for(tbb::blocked_range<int>(0, n), [&](const tbb::blocked_range<int>& r) {
        for (int j = r.begin(); j != r.end(); ++j) {
            for (int i = 0; i < null_dim; ++i) {
                NS_col[j][i] = NS[i][j];
            }
        }
        });

    return NS_col;
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
                throw std::invalid_argument("NullSpace's G-S Process Error: Isotropic vector encountered (v*v = 0 mod P)");
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
        throw std::invalid_argument("Ax=b calculation Error: Size is different");
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
        throw std::invalid_argument("Ax=b calculation Error: Size is different");
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
            throw std::invalid_argument("Matrix Diagonalize Error: Discriminant is a non-residue. Eigenvalues exist in F_p^2.");
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
            throw std::invalid_argument("Matrix Diagonalize Error: 2x2 matrix is defective (not diagonalizable)");
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
            throw std::invalid_argument("Matrix Diagonalize Error: 4x4 matrix is defective (not diagonalizable)");
        }
    }
}
inline void matrix_diagonalize_BF(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    int n = A.size();
    if (n == 0 || A.front().size() != n) {
        throw std::invalid_argument("Matrix diagonalization Error: Matrix is not square");
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
    throw std::invalid_argument("Matrix Diagonalize Error: Matrix is defective (not diagonalizable over F_p)");
}


















__global__ void build_S_and_SD_kernel(const long long* R, const long long* diag, long long* S, long long* SD, int n, long long mod) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n && j < n) {
        unsigned long long local_sum = 0;
        int limit = (i < j) ? i : j; // min(i, j)

        for (int k = 0; k <= limit; ++k) {
            // L[i][k]: 1 on diagonal, R[i][k] below diagonal
            long long Lik = (k == i) ? 1 : R[i * n + k];

            // U[k][j]: R[k][j] on and above diagonal (with 0-check on diagonal)
            long long Ukj = R[k * n + j];
            if (k == j && Ukj == 0) Ukj = 1;

            local_sum += (unsigned long long)Lik * Ukj;
            if (local_sum >> 61) local_sum %= mod; // Deferred modulo trick
        }

        long long S_ij = local_sum % mod;
        S[i * n + j] = S_ij;
        SD[i * n + j] = (S_ij * diag[j]) % mod;
    }
}

__global__ void matmul_kernel(const long long* A, const long long* B, long long* C, int n, long long mod) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n && j < n) {
        unsigned long long local_sum = 0;
        for (int k = 0; k < n; ++k) {
            local_sum += (unsigned long long)A[i * n + k] * B[k * n + j];
            if (local_sum >> 61) local_sum %= mod;
        }
        C[i * n + j] = local_sum % mod;
    }
}

inline vector<vector<long long>> generate_test_matrix_gpu(int N, long long MOD, bool sparse_eigenvalues = false, bool invertible = false) {
    // ---------------------------------------------------------
    // 1. Thread-Safe CPU RNG (Sequential)
    // ---------------------------------------------------------
    vector<long long> flat_R(N * N);
    vector<long long> diag(N);

    int pool_size = max(1, N / 100);
    vector<long long> eigen_pool(pool_size);
    if (sparse_eigenvalues) {
        for (int i = 0; i < pool_size; ++i) {
            long long val = get_rand(MOD);
            if (invertible && val == 0) val = 1;
            eigen_pool[i] = val;
        }
    }

    for (int i = 0; i < N * N; ++i) {
        flat_R[i] = get_rand(MOD);
    }

    for (int i = 0; i < N; ++i) {
        if (sparse_eigenvalues) {
            diag[i] = eigen_pool[get_rand(MOD) % pool_size];
        }
        else {
            long long val = get_rand(MOD);
            if (invertible && val == 0) val = 1;
            diag[i] = val;
        }
    }

    // ---------------------------------------------------------
    // 2. Async Allocations & Setup
    // ---------------------------------------------------------
    cudaStream_t stream;
    cudaStreamCreate(&stream); // Dedicated stream for TBB concurrency

    size_t bytes = N * N * sizeof(long long);
    long long* d_R, * d_diag, * d_S, * d_SD, * d_S_inv, * d_Final;
    long long* d_mults, * d_pivot_val; // [FIX]: Added d_pivot_val
    int* d_pivot_row, * d_info;
    int* h_info; // Pinned memory for error checking

    cudaMallocAsync(&d_R, bytes, stream);
    cudaMallocAsync(&d_diag, N * sizeof(long long), stream);
    cudaMallocAsync(&d_S, bytes, stream);
    cudaMallocAsync(&d_SD, bytes, stream);
    cudaMallocAsync(&d_S_inv, bytes, stream);
    cudaMallocAsync(&d_Final, bytes, stream);

    cudaMallocAsync(&d_mults, N * sizeof(long long), stream);
    cudaMallocAsync(&d_pivot_row, sizeof(int), stream);
    cudaMallocAsync(&d_pivot_val, sizeof(long long), stream); // [FIX]: Allocate d_pivot_val
    cudaMallocAsync(&d_info, sizeof(int), stream);

    cudaMallocHost(&h_info, sizeof(int));

    // Zero out the info flag on the device
    cudaMemsetAsync(d_info, 0, sizeof(int), stream);

    // Copy initial data to GPU asynchronously
    cudaMemcpyAsync(d_R, flat_R.data(), bytes, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_diag, diag.data(), N * sizeof(long long), cudaMemcpyHostToDevice, stream);

    // GPU Launch Configurations
    dim3 threads2D(16, 16);
    dim3 blocks2D((N + 15) / 16, (N + 15) / 16);
    int threads1D = 256;
    int blocks1D = (N + threads1D - 1) / threads1D;

    // ---------------------------------------------------------
    // 3. 100% Asynchronous GPU Pipeline (No CPU Syncs)
    // ---------------------------------------------------------

    // Build S and SD simultaneously from R
    build_S_and_SD_kernel << <blocks2D, threads2D, 0, stream >> > (d_R, d_diag, d_S, d_SD, N, MOD);

    // Invert d_S natively on the GPU 
    initIdentityKernel << <blocks2D, threads2D, 0, stream >> > (d_S_inv, N);

    for (int p = 0; p < N; ++p) {
        // [FIX]: Pass d_pivot_val so the scalar is recorded before normalization
        findPivotKernel << <1, 256, 0, stream >> > (
            d_S, N, N, p, nullptr, p, d_pivot_row, d_pivot_val, d_info);

        // Update singularity status on device
        updateDetInfoKernel << <1, 1, 0, stream >> > (
            d_pivot_row, nullptr, p, N, nullptr, d_info, MOD);

        swapRowsKernel<true> << <blocks1D, threads1D, 0, stream >> > (
            d_S, d_S_inv, N, N, N, p, nullptr, d_pivot_row, d_info);

        // [FIX]: Pass d_pivot_val to normalize without block-level race hazards
        normalizePivotRowKernel<true> << <blocks1D, threads1D, 0, stream >> > (
            d_S, d_S_inv, N, N, N, p, nullptr, p, d_pivot_row, MOD, d_pivot_val, d_info);

        computeMultipliersKernel << <blocks1D, threads1D, 0, stream >> > (
            d_S, N, N, p, nullptr, p, d_pivot_row, MOD, d_mults, true, d_info);

        eliminateRowsKernel<true> << <blocks2D, threads2D, 0, stream >> > (
            d_S, d_S_inv, N, N, N, p, nullptr, d_pivot_row, d_mults, MOD, d_info);
    }

    // Final Multiplication: A = SD * S_inv
    matmul_kernel << <blocks2D, threads2D, 0, stream >> > (d_SD, d_S_inv, d_Final, N, MOD);

    // ---------------------------------------------------------
    // 4. Single Synchronization Point
    // ---------------------------------------------------------
    vector<long long> flat_Final(N * N);

    cudaMemcpyAsync(h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync(flat_Final.data(), d_Final, bytes, cudaMemcpyDeviceToHost, stream);

    cudaStreamSynchronize(stream); // CPU waits ONLY here

    if (*h_info == 1) {
        throw std::invalid_argument("Error in Generation: Generated S is singular");
    }

    // ---------------------------------------------------------
    // 5. Cleanup
    // ---------------------------------------------------------
    cudaFreeAsync(d_R, stream);
    cudaFreeAsync(d_diag, stream);
    cudaFreeAsync(d_S, stream);
    cudaFreeAsync(d_SD, stream);
    cudaFreeAsync(d_S_inv, stream);
    cudaFreeAsync(d_Final, stream);
    cudaFreeAsync(d_mults, stream);
    cudaFreeAsync(d_pivot_row, stream);
    cudaFreeAsync(d_pivot_val, stream); // [FIX]: Free d_pivot_val
    cudaFreeAsync(d_info, stream);

    cudaFreeHost(h_info);
    cudaStreamDestroy(stream);

    // Convert to 2D vector for return
    vector<vector<long long>> Out(N, vector<long long>(N));
    for (int i = 0; i < N; ++i) {
        copy(flat_Final.begin() + i * N, flat_Final.begin() + (i + 1) * N, Out[i].begin());
    }

    return Out;
}


inline vector<vector<long long>> generate_test_matrix(int N, long long MOD, bool sparse_eigenvalues = false, bool invertible = false) {
    // 1. Thread-Safe RNG: Pre-generate all random numbers sequentially.
    // Memory allocation and a simple 2D sequential loop is incredibly fast 
    // and completely sidesteps any RNG race conditions.
    vector<vector<long long>> R(N, vector<long long>(N));
    vector<long long> diag(N);

    // Pre-generate a small pool for sparse (high-multiplicity) eigenvalues
    int pool_size = N / 100;
	if (pool_size < 1) pool_size = 1;
    vector<long long> eigen_pool(pool_size);
    if (sparse_eigenvalues) {
        for (int i = 0; i < pool_size; ++i) {
            long long val = get_rand(MOD);
            if (invertible && val == 0) val = 1; // Prevent 0 eigenvalue
            eigen_pool[i] = val;
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            R[i][j] = get_rand(MOD);
        }

        // Choose eigenvalue distribution based on the flag
        if (sparse_eigenvalues) {
            // Pick from the limited pool to guarantee repeating eigenvalues
            diag[i] = eigen_pool[get_rand(MOD) % pool_size];
        }
        else {
            // Dense: practically all distinct
            long long val = get_rand(MOD);
            if (invertible && val == 0) val = 1; // Prevent 0 eigenvalue
            diag[i] = val;
        }
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
            for (size_t j = 0; j < b.size() - 1; ++j) {
                long long sub = (factor * b[j]) % MOD;
                A[i - b.size() + 1 + j] = (A[i - b.size() + 1 + j] - sub + MOD) % MOD;
            }
            A[i] = 0; // Analytically guaranteed
        }
        trim(A);
        return A;
    }

    // O(N^2) Polynomial Long Division (Returns QUOTIENT A / B)
    vector<long long> div(vector<long long> A, const vector<long long>& B, long long MOD) {
        trim(A);
        vector<long long> b = B; trim(b);
        if (b.size() == 1 && b[0] == 0) return { 0 }; // Div by zero fallback
        if (A.size() < b.size()) return { 0 };

        vector<long long> Q(A.size() - b.size() + 1, 0);
        long long inv_lead = inverse(b.back());

        for (int i = A.size() - 1; i >= (int)b.size() - 1; --i) {
            if (A[i] == 0) continue;
            long long factor = (A[i] * inv_lead) % MOD;

            // Store in quotient
            Q[i - b.size() + 1] = factor;

            for (size_t j = 0; j < b.size() - 1; ++j) {
                long long sub = (factor * b[j]) % MOD;
                A[i - b.size() + 1 + j] = (A[i - b.size() + 1 + j] - sub + MOD) % MOD;
            }
            A[i] = 0; // Analytically guaranteed
        }
        trim(Q);
        return Q;
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

        vector<long long> next_v(m, 0);
        for (int i = 0; i <= min(m - 1, j + 1); ++i) {
            for (int k = max(0, i - 1); k <= min(m - 1, j); ++k) {
                next_v[i] = (next_v[i] + H[start + i][start + k] * curr[k]) % MOD;
            }
        }
        curr = next_v;
    }

    // Solve K * c = -H^m * e_1
    vector<long long> target(m, 0);
    for (int i = 0; i < m; ++i) target[i] = (MOD - curr[i]) % MOD;

    vector<long long> c(m, 0);
    for (int i = m - 1; i >= 0; --i) {
        unsigned long long sum = target[i];
        for (int j = i + 1; j < m; ++j) {
            // (MOD - K) handles the implicit subtraction: target - (K * c)
            sum += (unsigned long long)(MOD - K[i][j]) * c[j];
            if (sum >> 61) sum %= MOD; // Defer modulo safely
        }
        sum %= MOD;
        c[i] = (sum * inverse(K[i][i])) % MOD;
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
            total_poly = PolyMath::mul(total_poly, block_poly, MOD);
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

        // Due to Hessenberg structure, non-zeros only exist down to row = col + 1
        int limit = min(n, col + 2);

        while (sel < limit && M[sel][col] == 0) sel++;
        if (sel == limit) continue;

        if (sel != row) swap(M[row], M[sel]);
        pivot_col_to_row[col] = row;

        long long inv = inverse(M[row][col]);
        for (int j = col; j < n; ++j) M[row][j] = (M[row][j] * inv) % MOD;

        // OPTIMIZATION: Bounded elimination. We never have to eliminate rows below col + 1.
        for (int i = row + 1; i < limit; ++i) {
            if (M[i][col] != 0) {
                long long factor = M[i][col];
                for (int j = col + 1; j < n; ++j) {
                    long long sub = (factor * M[row][j]) % MOD;
                    M[i][j] = M[i][j] - sub + MOD;
					if (M[i][j] >= MOD) M[i][j] -= MOD;
                }
                M[i][col] = 0;
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
    vector<long long> a = { get_rand(MOD), 1 }; // a(x) = x + c

    // Compute d(x) = gcd(a(x)^((p-1)/2) - 1, poly(x))
    vector<long long> a_pow = PolyMath::pow_mod(a, (MOD - 1) / 2, poly, MOD);
    a_pow[0] = (a_pow[0] - 1 + MOD) % MOD; // Subtract 1

    vector<long long> d = PolyMath::gcd(a_pow, poly, MOD);

    // If it successfully split the polynomial into non-trivial factors, recurse
    if (d.size() > 1 && d.size() < poly.size()) {
        cz_split(d, roots, MOD);

        // BUG FIX: Use exact division (div) instead of modulo (mod)
        cz_split(PolyMath::div(poly, d, MOD), roots, MOD);
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

// 5. Updated Main Wrapper (Now with Lock-Free TBB Parallelization)
void matrix_diagonalize_krylov(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, long long MOD) {
    int n = A.size();
    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    // Phase 1: Reduce to Hessenberg (O(n^3))
    HessenbergResult hess = reduce_to_hessenberg_with_V(A, MOD);

    // Phase 2: Get characteristic polynomial from reduced H (O(n^3))
    vector<long long> poly = get_char_poly_krylov(hess.H, MOD);

    // Phase 3: Find eigenvalues
    vector<long long> roots = find_roots_Fp(poly, MOD);

    // Phase 4: Find Eigenvectors of H and map back to A in PARALLEL
    std::atomic<int> col_idx{ 0 };

    tbb::parallel_for(tbb::blocked_range<size_t>(0, roots.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i != r.end(); ++i) {
                // Early exit: Stop if the matrix is fully populated
                if (col_idx.load(std::memory_order_relaxed) >= n) return;

                long long lambda = roots[i];

                // O(n^2) nullspace extraction
                vector<vector<long long>> H_basis = get_hessenberg_nullspace(hess.H, lambda, MOD);

                for (const auto& y : H_basis) {
                    // 1. Atomically claim our column index BEFORE doing the math
                    // This entirely eliminates the need for std::mutex.
                    int current_col = col_idx.fetch_add(1, std::memory_order_relaxed);
                    
                    if (current_col >= n) return; // Over-allocated, stop work

                    // 2. Do the O(n^2) back-substitution directly into the matrix S
                    for (int r_idx = 0; r_idx < n; ++r_idx) {
                        
                        // 3. Use 128-bit integer to avoid inner-loop modulo branching 
                        // (Requires 64-bit GCC/Clang). If on 32-bit MSVC, revert to your original sum >> 61 logic.
                        long long sum = 0; 
                        
                        for (int c_idx = 0; c_idx < n; ++c_idx) {
                            if (y[c_idx]) {
                                sum += (long long)hess.V[r_idx][c_idx] * y[c_idx];
                                if (sum >> 61) sum %= MOD;
                            }
                        }
                        
                        // Modulo happens exactly once per row instead of randomly inside the loop
                        S[r_idx][current_col] = (long long)(sum % MOD);
                    }
                    
                    // Safely write the eigenvalue to the diagonal
                    D[current_col][current_col] = lambda;
                }
            }
        });
}

//// 5. Updated Main Wrapper (Sequential Version)
//void matrix_diagonalize_krylov(const vector<vector<long long>>& A, vector<vector<long long>>& S, vector<vector<long long>>& D, long long MOD) {
//    int n = A.size();
//    S.assign(n, vector<long long>(n, 0));
//    D.assign(n, vector<long long>(n, 0));
//
//    // Phase 1: Reduce to Hessenberg (O(n^3))
//    HessenbergResult hess = reduce_to_hessenberg_with_V(A, MOD);
//
//    // Phase 2: Get characteristic polynomial from reduced H (O(n^3))
//    vector<long long> poly = get_char_poly_krylov(hess.H, MOD);
//
//    // Phase 3: Find eigenvalues
//    vector<long long> roots = find_roots_Fp(poly, MOD);
//
//    // Phase 4: Find Eigenvectors of H and map back to A sequentially
//    int col_idx = 0;
//
//    for (size_t i = 0; i < roots.size(); ++i) {
//        // Early exit: Stop if the matrix is fully populated
//        if (col_idx >= n) break;
//
//        long long lambda = roots[i];
//
//        // O(n^2) nullspace extraction
//        vector<vector<long long>> H_basis = get_hessenberg_nullspace(hess.H, lambda, MOD);
//
//        for (const auto& y : H_basis) {
//            // 1. Claim our column index
//            int current_col = col_idx++;
//
//            if (current_col >= n) break; // Over-allocated, stop work
//
//            // 2. Do the O(n^2) back-substitution directly into the matrix S
//            for (int r_idx = 0; r_idx < n; ++r_idx) {
//
//                // 3. Use 128-bit integer to avoid inner-loop modulo branching 
//                // (Requires 64-bit GCC/Clang). If on 32-bit MSVC, revert to your original sum >> 61 logic.
//                long long sum = 0;
//
//                for (int c_idx = 0; c_idx < n; ++c_idx) {
//                    if (y[c_idx]) {
//                        sum += (long long)hess.V[r_idx][c_idx] * y[c_idx];
//                        if (sum >> 61) sum %= MOD;
//                    }
//                }
//
//                // Modulo happens exactly once per row instead of randomly inside the loop
//                S[r_idx][current_col] = (long long)(sum % MOD);
//            }
//
//            // Safely write the eigenvalue to the diagonal
//            D[current_col][current_col] = lambda;
//        }
//
//        // Break outer loop if matrix is filled during inner loop iterations
//        if (col_idx >= n) break;
//    }
//}









inline void matrix_diagonalize_henry(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
    int n = (int)A.size();
    S.assign(n, vector<long long>(n, 0));
    D.assign(n, vector<long long>(n, 0));

    vector<vector<long long>> AP_1 = matrix_power_gpu(A, MOD - 1);
    vector<vector<long long>> ZN1 = Null_Space_gpu(AP_1 - I_n(n), Orth);
    int eigvec_count = ZN1.empty() ? 0 : (int)ZN1[0].size();
    if (eigvec_count > 0)
        for (int j = 0; j < n; ++j)
            copy(ZN1[j].begin(), ZN1[j].end(), S[j].begin());
    vector<vector<long long>> ZN2 = Null_Space_gpu(AP_1, Orth);
    if (!ZN2.empty())
        for (int j = 0; j < n; ++j)
            copy(ZN2[j].begin(), ZN2[j].end(), S[j].begin() + eigvec_count);

    vector<vector<long long>> New_A(eigvec_count, vector<long long>(eigvec_count, 0));
    if (eigvec_count > 0 && eigvec_count < n) {
        vector<vector<long long>> S_inv = matrix_inverse_gpu(S);
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

    int mat_i = 0;

    for (int pi = 0; pi < MOD_decompose.size(); ++pi) {
        int start_mat_i = mat_i;
        int mati_upperbound = (int)M.size();

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
                    int current_col = 0;

                    vector<vector<long long>> PM = matrix_power_gpu(M[m_idx], powC);
                    vector<vector<long long>> St(m_size, vector<long long>(m_size, 0));
                    vector<int> eigspace_dim;

                    long long seed = seeds[FE[m_idx]] * inverse(MOD_decompose[pi]) % MOD;
                    long long seed2 = power(primitive, seed);

                    struct RootResult {
                        long long candidate;
                        vector<vector<long long>> ZN;
                    };
                    vector<RootResult> valid_roots;

                    std::mutex roots_mtx;
                    atomic<int> local_eigvec_count{ 0 };
                    if (roots_size <= 5) {
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

                                    vector<vector<long long>> ZN = query.size() > 30 ? Null_Space_gpu(query, Orth) : Null_Space(query, Orth);

                                    if (!ZN.empty()) {
                                        std::lock_guard<std::mutex> lock(roots_mtx);
                                        if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;
                                        valid_roots.push_back({ candidate, ZN });
                                        local_eigvec_count.fetch_add(ZN[0].size(), std::memory_order_relaxed);
                                    }
                                }
                            });
                    }
                    else {
                        HessenbergResult hess_PM = reduce_to_hessenberg_with_V(PM, MOD);
                        vector<long long> fp_roots = find_roots_Fp(get_char_poly_krylov(hess_PM.H, MOD), MOD);

                        tbb::parallel_for(tbb::blocked_range<size_t>(0, fp_roots.size()),
                            [&](const tbb::blocked_range<size_t>& inner_r) {
                                for (size_t i = inner_r.begin(); i != inner_r.end(); ++i) {
                                    // Early exit: Stop if other threads already found the full basis
                                    if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;

                                    long long candidate = fp_roots[i];

                                    // O(n^2) nullspace extraction on the upper Hessenberg matrix
                                    vector<vector<long long>> H_basis = get_hessenberg_nullspace(hess_PM.H, candidate, MOD);

                                    if (!H_basis.empty()) {
                                        int k = H_basis.size();
                                        vector<vector<long long>> ZN = hess_PM.V * matrix_transpose_tiled(H_basis);
                                        std::lock_guard<std::mutex> lock(roots_mtx);
                                        if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;
                                        valid_roots.push_back({ candidate, std::move(ZN) });
                                        local_eigvec_count.fetch_add(k, std::memory_order_relaxed);
                                    }
                                }
                            });
                    }

                    for (auto& root : valid_roots) {
                        if (current_col >= m_size) break;

                        local_result.new_FEs.push_back(root.candidate);
                        eigspace_dim.push_back((int)root.ZN[0].size());

                        for (int j = 0; j < root.ZN.size(); ++j)
                            copy(root.ZN[j].begin(), root.ZN[j].end(), St[j].begin() + current_col);
                        current_col += root.ZN[0].size();
                    }

                    vector<vector<long long>> mt = St.size() > 30 ? matrix_inverse_gpu(St) * M[m_idx] * St : matrix_inverse(St) * M[m_idx] * St;
                    for (int i = 0; i < St.size(); ++i)
                        copy(St[i].begin(), St[i].end(), ST[i + local_stp].begin() + local_stp);
                    matrix_chop(local_result.new_matrices, mt, eigspace_dim);
                }
            });

        mat_i = mati_upperbound;
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

inline void matrix_diagonalize_henry_naive(vector<vector<long long>> A, vector<vector<long long>>& S, vector<vector<long long>>& D, bool Orth) {
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

                    if (roots_size > 5) {
                        vector<vector<long long>> D, S;
                        matrix_diagonalize_krylov(M[m_idx], S, D, MOD);

                        for (int r = 0; r < N; ++r) {
                            local_result.new_matrices.push_back({ {D[r][r]} });
                            local_result.new_FEs.push_back(D[r][r]);
                            for (int c = 0; c < N; ++c) {
                                ST[local_stp + r][local_stp + c] = S[r][c];
                            }
                        }
                        continue;
                    }


                    vector<vector<long long>> PM = matrix_power(M[m_idx], powC);
                    vector<vector<long long>> St(m_size, vector<long long>(m_size, 0));
                    vector<int> eigspace_dim;
                    int current_col = 0;
                    //int local_eigvec_count = 0;

                    long long seed = seeds[FE[m_idx]] * inverse(MOD_decompose[pi]) % MOD;
                    long long seed2 = power(primitive, seed);

                    struct RootResult {
                        long long candidate;
                        vector<vector<long long>> ZN;
                    };
                    vector<RootResult> valid_roots;

                    std::mutex roots_mtx;
                    atomic<int> local_eigvec_count{ 0 };

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
                                    std::lock_guard<std::mutex> lock(roots_mtx);
                                    if (local_eigvec_count.load(std::memory_order_relaxed) >= m_size) return;
                                    valid_roots.push_back({ candidate, ZN });
                                    local_eigvec_count.fetch_add(ZN[0].size(), std::memory_order_relaxed);
                                }
                            }
                        });

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
















// --- Experiment B: Full Traditional vs Full Proposed Algorithm ---
inline void EX_B() {
    // 1. Create a dynamic filename based on the MOD value
    std::string filename = "experiment_B_results_" + std::to_string(MOD) + ".txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        printf("Error: Could not open %s for writing.\n", filename.c_str());
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
    //vector<int> N_values = { 512, 1024, 2048, 4096, 8192 };
    vector<int> N_values;
    for (int i = 1; i < 10; ++i)
        N_values.push_back(i * 300);

    print_and_write("=========================================================================================================\n");
    // Dynamically insert the MOD value using std::to_string
    print_and_write("Experiment B: Full Krylov Diagonalization vs Proposed Algorithm (" + std::to_string(MOD) + ")\n");
    print_and_write("=========================================================================================================\n");
    print_and_write("N\tTrials\tTrad Time (s)\tTrad Exp (x)\tProp Time (s)\tProp Exp (x)\tSpeedup\n");
    print_and_write("---------------------------------------------------------------------------------------------------------\n");

    int prev_N = 0;
    double prev_avg_krylov = 0.0;
    double prev_avg_henry = 0.0;

    for (int N : N_values) {
        int num_trials = (N <= 2000) ? 10 : 5;

        double total_time_krylov = 0.0;
        double total_time_henry = 0.0;

        // Ensure both algorithms have fresh output matrices
        std::vector<std::vector<long long>> S_krylov, D_krylov;
        std::vector<std::vector<long long>> S_henry, D_henry;

        for (int trial = 0; trial < num_trials; ++trial) {
            printf("G");
            std::vector<std::vector<long long>> A = generate_test_matrix_gpu(N, MOD, false, false);

            // --- Time the Full Traditional Baseline ---
            printf("K");
            auto start_krylov = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_krylov(A, S_krylov, D_krylov, MOD);
            auto end_krylov = std::chrono::high_resolution_clock::now();

            if (A * S_krylov != S_krylov * D_krylov) {
                printf("Krylov Diagonalization Error: A*S != S*D\n");
                exit(1);
            }
			if (matrix_rank(S_krylov) != N) {
				printf("Krylov Diagonalization Error: S is not full rank\n");
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
            if (matrix_rank(S_henry) != N) {
                printf("Henry Diagonalization Error: S is not full rank\n");
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

    // 4. Close the file and notify the user using the dynamic filename
    outfile.close();
    printf("Results successfully saved to '%s'\n", filename.c_str());
}

// --- Experiment C: Parallel Scaling Evaluation (O(I(n)) Proof) ---
inline void EX_C() {
    // 1. Create a dynamic filename based on the MOD value
    std::string filename = "experiment_C_results_" + std::to_string(MOD) + ".txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        printf("Error: Could not open %s for writing.\n", filename.c_str());
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

    // This array can now hold ANY arbitrary sequence of N values
    //vector<int> N_values = { 32, 64, 128, 256, 512, 1024, 2048, 4096 };
    vector<int> N_values;
    for (int i = 0; i < 10; ++i)
        N_values.push_back((i + 1) * 300);

    print_and_write("=========================================================================================================\n");
    // Dynamically insert the MOD value using std::to_string
    print_and_write("Experiment C: Parallel Diagonalization vs Parallel Inversion Scaling (" + std::to_string(MOD) + ")\n");
    print_and_write("=========================================================================================================\n");
    print_and_write("N\tTrials\tInversion (s)\tInv Exp (x)\tDiag Total (s)\tDiag Exp (x)\tRatio (Diag/Inv)\n");
    print_and_write("---------------------------------------------------------------------------------------------------------\n");

    double prev_inv_time = -1.0;
    double prev_diag_time = -1.0;
    int prev_N = -1; // Added to track the previous matrix dimension

    for (int N : N_values) {
        // Fewer trials for very large matrices
        int num_trials = (N <= 2000) ? 10 : 5;

        double total_time_inv = 0.0;
        double total_time_diag = 0.0;

        vector<vector<long long>> S_out, D_out;

        for (int trial = 0; trial < num_trials; ++trial) {
            vector<vector<long long>> A = generate_test_matrix_gpu(N, MOD, false, true);

            // --- Time Parallel Matrix Inversion I(n) ---
            auto start_inv = std::chrono::high_resolution_clock::now();
            vector<vector<long long>> A_inv = matrix_inverse(A);
            auto end_inv = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> duration_inv = end_inv - start_inv;
            total_time_inv += duration_inv.count();

            // --- Time Parallel Proposed Diagonalization ---
            auto start_diag = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_henry(A, S_out, D_out, false);
            auto end_diag = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> duration_diag = end_diag - start_diag;
            total_time_diag += duration_diag.count();
        }

        double avg_inv = total_time_inv / num_trials;
        double avg_diag = total_time_diag / num_trials;

        // Calculate empirical exponents dynamically using the general formula
        double exp_inv = 0.0;
        double exp_diag = 0.0;

        if (prev_inv_time > 0.0 && prev_N > 0) {
            double n_ratio = (double)N / prev_N;
            // x = ln(T_new / T_old) / ln(N_new / N_old)
            exp_inv = log(avg_inv / prev_inv_time) / log(n_ratio);
            exp_diag = log(avg_diag / prev_diag_time) / log(n_ratio);
        }

        double ratio = avg_diag / avg_inv;

        // 2. Format the numbers safely into a buffer
        char buffer[256];
        if (prev_inv_time < 0.0) {
            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\t-\t\t%lf\t-\t\t%.2fx\n", N, num_trials, avg_inv, avg_diag, ratio);
        }
        else {
            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\t%.3f\t\t%lf\t%.3f\t\t%.2fx\n", N, num_trials, avg_inv, exp_inv, avg_diag, exp_diag, ratio);
        }

        // 3. Send formatted buffer to both console and file
        print_and_write(buffer);

        // Store current metrics for the next iteration's exponent calculation
        prev_inv_time = avg_inv;
        prev_diag_time = avg_diag;
        prev_N = N;
    }
    print_and_write("=========================================================================================================\n");

    // 4. Close the file and notify the user using the dynamic filename
    outfile.close();
    printf("Results successfully saved to '%s'\n", filename.c_str());
}

// --- Experiment D: Naive Henry vs Optimized Henry Algorithm ---
inline void EX_D() {
    // 1. Open the text file for writing
    std::ofstream outfile("experiment_D_results.txt");
    if (!outfile.is_open()) {
        printf("Error: Could not open experiment_D_results.txt for writing.\n");
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

    vector<int> N_values;
    for (int i = 1; i < 10; ++i)
        N_values.push_back(i * 500);

    print_and_write("=========================================================================================================\n");
    print_and_write("Experiment D: Naive Henry vs Optimized Henry Algorithm (Fixed p)\n");
    print_and_write("=========================================================================================================\n");
    print_and_write("N\tTrials\tNaive Time(s)\tNaive Exp(x)\tOpt Time (s)\tOpt Exp (x)\tSpeedup\n");
    print_and_write("---------------------------------------------------------------------------------------------------------\n");

    int prev_N = 0;
    double prev_avg_naive = 0.0;
    double prev_avg_opt = 0.0;

    for (int N : N_values) {
        int num_trials = (N <= 256) ? 10 : (N <= 512) ? 5 : 2;

        double total_time_naive = 0.0;
        double total_time_opt = 0.0;

        // Ensure both algorithms have fresh output matrices
        std::vector<std::vector<long long>> S_naive, D_naive;
        std::vector<std::vector<long long>> S_opt, D_opt;

        for (int trial = 0; trial < num_trials; ++trial) {
            printf("G");
            std::vector<std::vector<long long>> A = generate_test_matrix(N, MOD, true, false);

            // --- Time the Naive Henry Baseline ---
            printf("N");
            auto start_naive = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_henry_naive(A, S_naive, D_naive, false);
            auto end_naive = std::chrono::high_resolution_clock::now();

            if (A * S_naive != S_naive * D_naive) {
                printf("Naive Henry Error: A*S != S*D\n");
                exit(1);
            }
            if (matrix_rank(S_naive) != N) {
                printf("Naive Henry Error: S is not full rank\n");
                exit(1);
            }

            std::chrono::duration<double> duration_naive = end_naive - start_naive;
            total_time_naive += duration_naive.count();

            // --- Time the Optimized Henry Algorithm ---
            printf("O");
            auto start_opt = std::chrono::high_resolution_clock::now();
            matrix_diagonalize_henry(A, S_opt, D_opt, false);
            auto end_opt = std::chrono::high_resolution_clock::now();

            if (A * S_opt != S_opt * D_opt) {
                printf("Optimized Henry Error: A*S != S*D\n");
                exit(1);
            }
            if (matrix_rank(S_opt) != N) {
                printf("Optimized Henry Error: S is not full rank\n");
                exit(1);
            }

            std::chrono::duration<double> duration_opt = end_opt - start_opt;
            total_time_opt += duration_opt.count();
        }
        printf("\n");

        double avg_naive = total_time_naive / num_trials;
        double avg_opt = total_time_opt / num_trials;
        double speedup = avg_naive / avg_opt;

        // 2. Format the numbers safely into a buffer to retain your exact formatting
        char buffer[256];

        if (prev_N > 0 && prev_avg_naive > 0.0 && prev_avg_opt > 0.0) {
            // Calculate exponents
            double exp_naive = std::log(avg_naive / prev_avg_naive) / std::log((double)N / prev_N);
            double exp_opt = std::log(avg_opt / prev_avg_opt) / std::log((double)N / prev_N);

            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%.2fx\n",
                N, num_trials, avg_naive, exp_naive, avg_opt, exp_opt, speedup);
        }
        else {
            // First run, no previous data
            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\tN/A\t\t%lf\tN/A\t\t%.2fx\n",
                N, num_trials, avg_naive, avg_opt, speedup);
        }

        // 3. Send formatted buffer to both console and file
        print_and_write(buffer);

        // Store current values for the next iteration
        prev_N = N;
        prev_avg_naive = avg_naive;
        prev_avg_opt = avg_opt;
    }

    print_and_write("=========================================================================================================\n");

    // 4. Close the file and notify the user
    outfile.close();
    printf("Results successfully saved to 'experiment_D_results.txt'\n");
}

// --- Experiment B: Full Traditional vs Full Proposed Algorithm vs Inverse ---
inline void EX_FULL() {
    // 1. Create a dynamic filename based on the MOD value
    std::string filename = "experiment_FULL_results_" + std::to_string(MOD) + ".txt";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        printf("Error: Could not open %s for writing.\n", filename.c_str());
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

    vector<int> N_values;
    for (int i = 1; i <= 10; ++i)
        N_values.push_back(i * 300);
    //vector<int> N_values = { 2100, 2400, 2700, 3000 };

    print_and_write("=======================================================================================================================================\n");
    // Dynamically insert the MOD value using std::to_string
    print_and_write("Experiment FULL: Proposed (Henry) vs Krylov Diagonalization vs Matrix Inverse (" + std::to_string(MOD) + ")\n");
    print_and_write("=======================================================================================================================================\n");
    print_and_write("N\tTrials\tHenry(s)\tHenry(x)\t|\tKrylov(s)\tKrylov(x)\tSpeedup\t|\tInverse(s)\tInverse(x)\tSpeedup\n");
    print_and_write("---------------------------------------------------------------------------------------------------------------------------------------\n");

    int prev_N = 0;
    double prev_avg_henry = 0.0;
    double prev_avg_krylov = 0.0;
    double prev_avg_inverse = 0.0;

    for (int N : N_values) {
        int num_trials = (N <= 2000) ? 10 : 5;

        double total_time_henry = 0.0;
        double total_time_krylov = 0.0;
        double total_time_inverse = 0.0;

        for (int trial = 0; trial < num_trials; ++trial) {
            std::vector<std::vector<long long>> A;

            // 1. Guard Matrix Generation Separately
            bool matrix_generated = false;
            while (!matrix_generated) {
                try {
                    printf("G");
                    A = generate_test_matrix_gpu(N, MOD, false, true);
                    matrix_generated = true;
                }
                catch (const std::exception& e) {
                    printf("\n[Matrix Gen Error] %s - Retrying matrix generation...\n", e.what());
                }
                catch (...) {
                    printf("\n[Matrix Gen Error] Unknown exception - Retrying matrix generation...\n");
                }
            }

            // 2. Guard Algorithms and retry with the SAME matrix if they fail
            bool trial_successful = false;
            while (!trial_successful) {
                try {
                    // --- Time the Proposed Algorithm (Henry) ---
                    printf("H");
                    std::vector<std::vector<long long>> S_henry, D_henry;
                    auto start_henry = std::chrono::high_resolution_clock::now();
                    matrix_diagonalize_henry(A, S_henry, D_henry, false);
                    auto end_henry = std::chrono::high_resolution_clock::now();

                    if (A * S_henry != S_henry * D_henry) {
                        throw std::runtime_error("Henry Diagonalization Error: A*S != S*D");
                    }
                    if (matrix_rank_gpu(S_henry) != N) {
                        throw std::runtime_error("Henry Diagonalization Error: S is not full rank");
                    }

                    std::chrono::duration<double> duration_henry = end_henry - start_henry;

                    // --- Time the Full Traditional Baseline (Krylov) ---
                    printf("K");
                    std::vector<std::vector<long long>> S_krylov, D_krylov;
                    auto start_krylov = std::chrono::high_resolution_clock::now();
                    matrix_diagonalize_krylov(A, S_krylov, D_krylov, MOD);
                    auto end_krylov = std::chrono::high_resolution_clock::now();

                    if (A * S_krylov != S_krylov * D_krylov) {
                        throw std::runtime_error("Krylov Diagonalization Error: A*S != S*D");
                    }
                    if (matrix_rank_gpu(S_krylov) != N) {
                        throw std::runtime_error("Krylov Diagonalization Error: S is not full rank");
                    }

                    std::chrono::duration<double> duration_krylov = end_krylov - start_krylov;

                    // --- Time Matrix Inverse ---
                    printf("I");
                    auto start_inverse = std::chrono::high_resolution_clock::now();
                    std::vector<std::vector<long long>> A_inv = matrix_inverse_gpu(A);
                    auto end_inverse = std::chrono::high_resolution_clock::now();

                    std::chrono::duration<double> duration_inverse = end_inverse - start_inverse;

                    // ONLY add to total times if all three algorithms finished without throwing exceptions
                    total_time_henry += duration_henry.count();
                    total_time_krylov += duration_krylov.count();
                    total_time_inverse += duration_inverse.count();

                    // Mark as successful to break the retry loop
                    trial_successful = true;
                }
                catch (const std::exception& e) {
                    // WARNING: If the error is deterministic (e.g. matrix is singular), this will infinite loop.
                    printf("\n[Error] %s - Retrying trial with SAME matrix...\n", e.what());
                }
                catch (...) {
                    printf("\n[Error] Unknown exception occurred - Retrying trial with SAME matrix...\n");
                }
            }
        }
        printf("\n");

        double avg_henry = total_time_henry / num_trials;
        double avg_krylov = total_time_krylov / num_trials;
        double avg_inverse = total_time_inverse / num_trials;

        // Speedups are calculated relative to Henry
        double speedup_krylov = avg_krylov / avg_henry;
        double speedup_inverse = avg_inverse / avg_henry;

        // 2. Format the numbers safely into a larger buffer
        char buffer[512];

        if (prev_N > 0 && prev_avg_henry > 0.0 && prev_avg_krylov > 0.0 && prev_avg_inverse > 0.0) {
            // Calculate exponents
            double exp_henry = std::log(avg_henry / prev_avg_henry) / std::log((double)N / prev_N);
            double exp_krylov = std::log(avg_krylov / prev_avg_krylov) / std::log((double)N / prev_N);
            double exp_inverse = std::log(avg_inverse / prev_avg_inverse) / std::log((double)N / prev_N);

            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\t%lf\t|\t%lf\t%lf\t%.2fx\t|\t%lf\t%lf\t%.2fx\n",
                N, num_trials, avg_henry, exp_henry, avg_krylov, exp_krylov, speedup_krylov, avg_inverse, exp_inverse, speedup_inverse);
        }
        else {
            // First run, no previous data
            snprintf(buffer, sizeof(buffer), "%d\t%d\t%lf\tN/A\t\t|\t%lf\tN/A\t\t%.2fx\t|\t%lf\tN/A\t\t%.2fx\n",
                N, num_trials, avg_henry, avg_krylov, speedup_krylov, avg_inverse, speedup_inverse);
        }

        // 3. Send formatted buffer to both console and file
        print_and_write(buffer);

        // Store current values for the next iteration
        prev_N = N;
        prev_avg_henry = avg_henry;
        prev_avg_krylov = avg_krylov;
        prev_avg_inverse = avg_inverse;
    }

    print_and_write("=======================================================================================================================================\n");

    // 4. Close the file and notify the user using the dynamic filename
    outfile.close();
    printf("Results successfully saved to '%s'\n", filename.c_str());
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
            M[i][j] = power(v[j][0], i);

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
    for (int i = 0; i < m && b.size() < t + 1; ++i) {
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
    //MOD = 100000007;          //2*491*101833
    //MOD = 100663291;          //2*3*3*3*5*7*13*16*241
    //MOD = 131071;             //2*3*5*17*257
    //MOD = 524287;             //2*3*3*3*7*19*73
    //MOD = 65537;              //2^16
    //MOD = 653659;             //2*3*108943
    //MOD = 257;			    //2*2*2*2*2*2*2*2
    //MOD = 101;                //2*2*5*5

    //Safe primes
    //MOD = 565127;
    //MOD = 1000000007;


    //Initiation();

    //EX_B();
    //EX_C();



    vector<long long> mod_values = {
        100663291,  // 2*3*3*3*5*7*13*16*241
        131071,     // 2*3*5*17*257
        524287,     // 2*3*3*3*7*19*73

        100000007,  // 2*491*101833
        653659,     // 2*3*108943
        565127,     // Safe prime
    };

    for (long long current_mod : mod_values) {
        MOD = current_mod;

        printf("\n======================================================================\n");
        printf("Starting sequence for MOD = %lld\n", MOD);
        printf("======================================================================\n\n");

        Initiation();

		EX_FULL();
    }





 //   MOD = 653659;
	//Initiation();

	//vector<vector<long long>> M = generate_test_matrix_gpu(3000, MOD);
 //   
 //   auto start_cpu = std::chrono::high_resolution_clock::now();
 //   generate_test_matrix(3000, MOD);
 //   auto end_cpu = std::chrono::high_resolution_clock::now();
 //   std::chrono::duration<double, std::milli> time_cpu = end_cpu - start_cpu;

 //   printf("CPU Time taken: %.3f ms\n\n", time_cpu.count());

 //   auto start_gpu = std::chrono::high_resolution_clock::now();
 //   generate_test_matrix_gpu(3000, MOD);
 //   auto end_gpu = std::chrono::high_resolution_clock::now();
 //   std::chrono::duration<double, std::milli> time_gpu = end_gpu - start_gpu;

 //   printf("GPU Time taken: %.3f ms\n", time_gpu.count());


 //   vector<vector<long long>> M = generate_test_matrix(5, MOD, false, true);
	//vector<vector<long long>> D = I_n(5);
	//D[0][0] = 3, D[1][1] = 2, D[2][2] = 4, D[3][3] = 1, D[4][4] = 0;
 //   M = M * D * matrix_inverse(M);
	//vector<vector<long long>> MP_1 = matrix_power(M, MOD - 1);
 //   matrix_print(MP_1);
 //   matrix_print(Null_Space(MP_1 - I_n(5), false));
 //   matrix_print(Null_Space(MP_1, false));
 //   matrix_print(Null_Space(M, false));



    return 0;
}
