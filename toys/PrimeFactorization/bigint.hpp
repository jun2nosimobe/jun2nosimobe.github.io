#pragma once

#include <string>
#include <cstdlib>
#include <gmp.h>
#include <numeric>

// =================================================================
// 【mini-gmp C++ Wrapper】
// 既存の cpp_int の振る舞いを再現し、数式演算子をそのまま使えるようにする
// =================================================================
class cpp_int {
private:
    mpz_t val;
public:
    cpp_int() { mpz_init(val); }
    cpp_int(int v) { mpz_init_set_si(val, v); }
    cpp_int(long long v) { 
        std::string s = std::to_string(v);
        mpz_init_set_str(val, s.c_str(), 10); 
    }
    
    // ★ 追加: uint64_t (unsigned long long) 用のコンストラクタ
    cpp_int(unsigned long long v) { 
        std::string s = std::to_string(v);
        mpz_init_set_str(val, s.c_str(), 10); 
    }

    cpp_int(const char* s) { mpz_init_set_str(val, s, 10); }
    cpp_int(const std::string& s) { mpz_init_set_str(val, s.c_str(), 10); }
    
    cpp_int(const cpp_int& other) { mpz_init_set(val, other.val); }
    cpp_int(cpp_int&& other) noexcept { 
        val[0] = other.val[0]; 
        mpz_init(other.val); 
    }
    
    ~cpp_int() { mpz_clear(val); }

    cpp_int& operator=(const cpp_int& other) {
        if (this != &other) mpz_set(val, other.val);
        return *this;
    }
    cpp_int& operator=(int v) { mpz_set_si(val, v); return *this; }

    std::string str() const {
        char* s = mpz_get_str(NULL, 10, val);
        std::string res(s);
        free(s);
        return res;
    }

    explicit operator int() const { return (int)mpz_get_si(val); }
    explicit operator long long() const { 
        char* s = mpz_get_str(NULL, 10, val);
        long long res = std::stoll(s);
        free(s);
        return res;
    }

    const mpz_t& get_mpz_t() const { return val; }
    mpz_t& get_mpz_t() { return val; }

    // ★追加: テンプレートアルゴリズム（mont_ladder等）から呼ばれる共通インターフェース
    size_t bit_length() const { return mpz_sizeinbase(val, 2); }
    bool test_bit(int i) const { return mpz_tstbit(val, i); }

    // 算術演算子 (cpp_int同士)
    cpp_int operator+(const cpp_int& b) const { cpp_int r; mpz_add(r.val, val, b.val); return r; }
    cpp_int operator-(const cpp_int& b) const { cpp_int r; mpz_sub(r.val, val, b.val); return r; }
    cpp_int operator*(const cpp_int& b) const { cpp_int r; mpz_mul(r.val, val, b.val); return r; }
    cpp_int operator/(const cpp_int& b) const { cpp_int r; mpz_tdiv_q(r.val, val, b.val); return r; }
    cpp_int operator%(const cpp_int& b) const { cpp_int r; mpz_tdiv_r(r.val, val, b.val); return r; }
    cpp_int operator-() const { cpp_int r; mpz_neg(r.val, val); return r; } 
    
    // 算術演算子 (intとの計算)
    cpp_int operator+(int b) const { return *this + cpp_int(b); }
    cpp_int operator-(int b) const { return *this - cpp_int(b); }
    cpp_int operator*(int b) const { return *this * cpp_int(b); }
    cpp_int operator/(int b) const { return *this / cpp_int(b); }
    cpp_int operator%(int b) const { return *this % cpp_int(b); }

    friend cpp_int operator+(int a, const cpp_int& b) { return cpp_int(a) + b; }
    friend cpp_int operator-(int a, const cpp_int& b) { return cpp_int(a) - b; }
    friend cpp_int operator*(int a, const cpp_int& b) { return cpp_int(a) * b; }
    friend cpp_int operator/(int a, const cpp_int& b) { return cpp_int(a) / b; }
    friend cpp_int operator%(int a, const cpp_int& b) { return cpp_int(a) % b; }

    cpp_int operator<<(int b) const { cpp_int r; mpz_mul_2exp(r.val, val, b); return r; }
    cpp_int operator>>(int b) const { cpp_int r; mpz_fdiv_q_2exp(r.val, val, b); return r; }
    cpp_int operator&(const cpp_int& b) const { cpp_int r; mpz_and(r.val, val, b.val); return r; }
    cpp_int operator|(const cpp_int& b) const { cpp_int r; mpz_ior(r.val, val, b.val); return r; }

    cpp_int& operator+=(const cpp_int& b) { mpz_add(val, val, b.val); return *this; }
    cpp_int& operator-=(const cpp_int& b) { mpz_sub(val, val, b.val); return *this; }
    cpp_int& operator*=(const cpp_int& b) { mpz_mul(val, val, b.val); return *this; }
    cpp_int& operator/=(const cpp_int& b) { mpz_tdiv_q(val, val, b.val); return *this; }
    cpp_int& operator%=(const cpp_int& b) { mpz_tdiv_r(val, val, b.val); return *this; }
    cpp_int& operator>>=(int b) { mpz_fdiv_q_2exp(val, val, b); return *this; }
    cpp_int& operator<<=(int b) { mpz_mul_2exp(val, val, b); return *this; }

    cpp_int& operator+=(int b) { return *this += cpp_int(b); }
    cpp_int& operator-=(int b) { return *this -= cpp_int(b); }
    cpp_int& operator*=(int b) { return *this *= cpp_int(b); }
    cpp_int& operator/=(int b) { return *this /= cpp_int(b); }
    cpp_int& operator%=(int b) { return *this %= cpp_int(b); }

    bool operator==(const cpp_int& b) const { return mpz_cmp(val, b.val) == 0; }
    bool operator!=(const cpp_int& b) const { return mpz_cmp(val, b.val) != 0; }
    bool operator<(const cpp_int& b) const { return mpz_cmp(val, b.val) < 0; }
    bool operator>(const cpp_int& b) const { return mpz_cmp(val, b.val) > 0; }
    bool operator<=(const cpp_int& b) const { return mpz_cmp(val, b.val) <= 0; }
    bool operator>=(const cpp_int& b) const { return mpz_cmp(val, b.val) >= 0; }

    bool operator==(int b) const { return mpz_cmp_si(val, b) == 0; }
    bool operator!=(int b) const { return mpz_cmp_si(val, b) != 0; }
    bool operator<(int b) const { return mpz_cmp_si(val, b) < 0; }
    bool operator>(int b) const { return mpz_cmp_si(val, b) > 0; }
};

// 安全なモジュロ演算（負の数に対応）
inline cpp_int mod(cpp_int a, cpp_int n) {
    a %= n;
    if (a < 0) a += n;
    return a;
}

// Boostの gcd をエミュレートするグローバル関数
inline cpp_int gcd(const cpp_int& a, const cpp_int& b) {
    cpp_int r;
    mpz_gcd(r.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    return r;
}

// =================================================================
// 【384ビット専用の固定長整数】
// =================================================================
struct uint384_t {
    uint64_t limbs[6] = {0};

    uint384_t() {}
    uint384_t(uint64_t v) { limbs[0] = v; }
    
    // GMPの機能を使って文字列から高速にビット列へ変換
    uint384_t(const std::string& s) {
        mpz_t z;
        mpz_init_set_str(z, s.c_str(), 10);
        size_t count;
        mpz_export(limbs, &count, -1, sizeof(uint64_t), 0, 0, z);
        for(size_t i = count; i < 6; ++i) limbs[i] = 0;
        mpz_clear(z);
    }
    
    // ビット列から文字列への変換もGMPに任せる
    std::string str() const {
        mpz_t z;
        mpz_init(z);
        mpz_import(z, 6, -1, sizeof(uint64_t), 0, 0, limbs);
        char* s = mpz_get_str(NULL, 10, z);
        std::string res(s);
        free(s);
        mpz_clear(z);
        return res;
    }

    // 比較演算子
    bool operator==(const uint384_t& b) const {
        for(int i = 0; i < 6; i++) if(limbs[i] != b.limbs[i]) return false;
        return true;
    }
    bool operator==(uint64_t b) const {
        if (limbs[0] != b) return false;
        for(int i = 1; i < 6; i++) if(limbs[i] != 0) return false;
        return true;
    }
    
    bool operator>=(const uint384_t& b) const {
        for(int i = 5; i >= 0; i--) {
            if(limbs[i] > b.limbs[i]) return true;
            if(limbs[i] < b.limbs[i]) return false;
        }
        return true;
    }
    bool operator<(const uint384_t& b) const { return !(*this >= b); }
    bool operator>(const uint384_t& b) const { return (*this >= b) && !(*this == b); }
    bool operator<=(const uint384_t& b) const { return !(*this > b); }
    bool operator!=(const uint384_t& b) const { return !(*this == b); }

    uint384_t& operator+=(const uint384_t& b) { *this = *this + b; return *this; }
    uint384_t& operator-=(const uint384_t& b) { *this = *this - b; return *this; }
    uint384_t& operator*=(const uint384_t& b);

    // テンプレート互換用メソッド
    size_t bit_length() const {
        for (int i = 5; i >= 0; --i) {
            if (limbs[i] != 0) return i * 64 + 64 - __builtin_clzll(limbs[i]);
        }
        return 0;
    }
    
    bool test_bit(int i) const {
        if (i < 0 || i >= 384) return false;
        return (limbs[i / 64] >> (i % 64)) & 1;
    }

    // 128ビット幅を用いたキャリー付き加算（ループ展開される）
    uint384_t operator+(const uint384_t& b) const {
        uint384_t res;
        uint64_t carry = 0;
        for (int i = 0; i < 6; ++i) {
            __uint128_t sum = (__uint128_t)limbs[i] + b.limbs[i] + carry;
            res.limbs[i] = (uint64_t)sum;
            carry = (uint64_t)(sum >> 64);
        }
        return res;
    }
    
    // ボロー付き減算
    uint384_t operator-(const uint384_t& b) const {
        uint384_t res;
        uint64_t borrow = 0;
        for (int i = 0; i < 6; ++i) {
            __uint128_t diff = (__uint128_t)limbs[i] - b.limbs[i] - borrow;
            res.limbs[i] = (uint64_t)diff;
            borrow = (uint64_t)(diff >> 127); // 負なら1
        }
        return res;
    }
};

// =================================================================
// ★ 超高速連携: メモリを直接コピーして型変換する
// =================================================================
inline cpp_int to_gmp(const uint384_t& a) {
    cpp_int res;
    mpz_import(res.get_mpz_t(), 6, -1, sizeof(uint64_t), 0, 0, a.limbs);
    return res;
}

inline uint384_t from_gmp(const cpp_int& a) {
    uint384_t res;
    size_t count;
    uint64_t temp[16] = {0}; // ★修正: 最大1024ビットまで受け止められる安全なバッファ
    mpz_export(temp, &count, -1, sizeof(uint64_t), 0, 0, a.get_mpz_t());
    for(size_t i = 0; i < 6 && i < count; ++i) res.limbs[i] = temp[i];
    return res;
}

// =================================================================
// ★ 高速・安全なモジュロ演算モジュール (アンダーフロー・オーバーフロー完全対応)
// =================================================================
template <typename IntType>
inline IntType add_mod(const IntType& a, const IntType& b, const IntType& n) {
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        return mod(a + b, n);
    } else {
        uint384_t res = a + b;
        if (res < a || res >= n) res = res - n;
        return res;
    }
}

template <typename IntType>
inline IntType sub_mod(const IntType& a, const IntType& b, const IntType& n) {
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        return mod(a - b, n);
    } else {
        if (a >= b) return a - b;
        return n - (b - a); // 符号なし整数の安全な引き算
    }
}

template <typename IntType>
inline IntType mul_mod(const IntType& a, const IntType& b, const IntType& n) {
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        return mod(a * b, n);
    } else {
        // 掛け算の膨張(768bit)をGMP内で安全に還元してから戻す
        return from_gmp(mod(to_gmp(a) * to_gmp(b), to_gmp(n)));
    }
}

// --- 残りのブリッジ演算子 ---
inline uint384_t operator*(const uint384_t& a, const uint384_t& b) {
    return from_gmp(to_gmp(a) * to_gmp(b));
}
inline uint384_t operator%(const uint384_t& a, const uint384_t& b) {
    return from_gmp(to_gmp(a) % to_gmp(b));
}
inline uint384_t operator/(const uint384_t& a, const uint384_t& b) {
    return from_gmp(to_gmp(a) / to_gmp(b));
}
inline uint384_t mod(uint384_t a, uint384_t n) {
    return from_gmp(mod(to_gmp(a), to_gmp(n)));
}
inline uint384_t gcd(const uint384_t& a, const uint384_t& b) {
    return from_gmp(gcd(to_gmp(a), to_gmp(b)));
}
inline uint384_t& uint384_t::operator*=(const uint384_t& b) {
    *this = *this * b;
    return *this;
}