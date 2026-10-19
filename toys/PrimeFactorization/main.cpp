#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <string>
#include <cstdlib>

// 分割したモジュール群を読み込む
#include "bigint.hpp"
#include "montgomery.hpp"
#include "polynomial.hpp"
#include "ecm_core.hpp"

using namespace emscripten;

// =================================================================
// ユーティリティ関数
// =================================================================

// 素数判定用の繰り返し二乗法
cpp_int powm(cpp_int base, cpp_int exp, cpp_int mod_val) {
    cpp_int res = 1;
    base = mod(base, mod_val);
    while (exp > 0) {
        if (exp % 2 == 1) res = mod(res * base, mod_val);
        base = mod(base * base, mod_val);
        exp /= 2;
    }
    return res;
}

// 巨大な多倍長整数用の乱数生成関数
cpp_int generate_random_base(cpp_int limit) {
    if (limit <= 0) return 0;
    cpp_int res = 0;
    cpp_int temp = limit;
    int chunks = 0;
    while(temp > 0) { temp >>= 30; chunks++; }
    for(int i = 0; i < chunks; i++) res = (res << 30) | (std::rand() & 0x3FFFFFFF);
    return res % limit;
}

// 本格版ミラー・ラビン素数判定
bool is_prime_cpp(const std::string& n_str) {
    cpp_int n(n_str);
    if (n < 2) return false;
    int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
    for (int p : small_primes) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    cpp_int d = n - 1;
    int s = 0;
    while (d % 2 == 0) { d /= 2; s++; }
    
    int k = 20; 
    for (int i = 0; i < k; i++) {
        cpp_int limit = n - 3;
        cpp_int a = 2 + generate_random_base(limit);
        cpp_int x = powm(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int r = 1; r < s; r++) {
            x = mod(x * x, n);
            if (x == n - 1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

// ポラード・ロー法
std::string pollard_rho(const std::string& n_str, int max_steps) {
    cpp_int n(n_str);
    if (n <= 1) return "";
    if (n % 2 == 0) return "2";
    
    cpp_int x = 2;
    cpp_int y = 2;
    cpp_int d = 1;
    cpp_int c = cpp_int(std::rand()) % (n - 1) + 1;

    int step = 0;
    while (d == 1 && step < max_steps) {
        x = mod(x * x + c, n);
        y = mod(y * y + c, n);
        y = mod(y * y + c, n);
        
        cpp_int diff = x > y ? x - y : y - x;
        d = gcd(diff, n);
        
        step++;
    }

    if (d > 1 && d < n) {
        return d.str(); 
    }
    return ""; 
}

// =================================================================
// ディスパッチャ (JavaScriptからのエントリーポイント)
// =================================================================
val run_ecm_batch(std::string n_str, int B1, int seed, int num_curves) {
    cpp_int N_test(n_str);
    
    // Nのビット数を判定
    size_t bits = N_test.bit_length();
    
    if (bits <= 384) {
        // ★ 将来 uint384_t の実装が完了したら以下に切り替えます
        return run_ecm_core<uint384_t, MontgomeryContext384>(n_str, B1, seed, num_curves);
        
        // 現在は従来通りの GMP エンジンへルーティング
        //return run_ecm_core<cpp_int, MontgomeryContextGMP>(n_str, B1, seed, num_curves);
    } else {
        // 384ビットを超える場合は従来の汎用 GMP エンジンへルーティング
        return run_ecm_core<cpp_int, MontgomeryContextGMP>(n_str, B1, seed, num_curves);
    }
}

// =================================================================
// Wasm バインディング
// =================================================================
EMSCRIPTEN_BINDINGS(ecm_module) {
    function("is_prime_cpp", &is_prime_cpp);
    function("run_ecm_batch", &run_ecm_batch);
    function("pollard_rho", &pollard_rho);
}