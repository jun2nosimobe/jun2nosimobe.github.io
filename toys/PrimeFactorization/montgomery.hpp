#pragma once

#include "bigint.hpp"
#include <vector>

// --- テンプレート化: モントゴメリ曲線用のX, Z座標のみの構造体 ---
template <typename IntType>
struct PointXZ {
    IntType X, Z;
};

// =================================================================
// 【384ビット専用 CIOS モントゴメリ乗算エンジン】
// =================================================================
struct MontgomeryContext384 {
    uint384_t n;
    uint64_t n_prime; 
    uint384_t r2; // ★ 追加: R^2 mod N を保持するキャッシュ

    MontgomeryContext384(const uint384_t& modulus) {
        n = modulus;
        uint64_t inv = n.limbs[0]; 
        
        // ★ 修正: 反復回数が足りず上位16ビットがゴミになっていた致命的バグを修正
        // 確実な64ビット精度（最大96ビット相当）を得るために 6回ループさせる
        for (int i = 0; i < 6; ++i) {
            inv *= 2 - n.limbs[0] * inv;
        }
        n_prime = -inv;

        // R = 2^384 なので R^2 = 2^768
        cpp_int N_gmp = to_gmp(n);
        cpp_int R2 = (cpp_int(1) << 768) % N_gmp; 
        r2 = from_gmp(R2);
    }

    uint384_t to_mont(const uint384_t& a) const {
        cpp_int A(a.str());
        cpp_int N(n.str());
        cpp_int R = (cpp_int(1) << 384);
        cpp_int res = (A * R) % N;
        return uint384_t(res.str());
    }

    uint384_t from_mont(const uint384_t& a) const {
        return mont_mul(a, uint384_t((uint64_t)1));
    }

    // ★ 究極の CIOS (Coarsely Integrated Operand Scanning) 乗算 (完全修正版)
    uint384_t mont_mul(const uint384_t& a, const uint384_t& b) const {
        uint384_t res;
        uint64_t t[8] = {0}; // ★ 修正1: 桁溢れを完全に受け止めるため 8 要素に拡張
        
        for (int i = 0; i < 6; ++i) {
            uint64_t c = 0;
            // 乗算ステップ
            for (int j = 0; j < 6; ++j) {
                __uint128_t prod = (__uint128_t)a.limbs[i] * b.limbs[j] + t[j] + c;
                t[j] = (uint64_t)prod;
                c = (uint64_t)(prod >> 64);
            }
            
            // ★ 修正2: t[6] += c; によるサイレント・オーバーフローを防止
            __uint128_t sum6 = (__uint128_t)t[6] + c;
            t[6] = (uint64_t)sum6;
            t[7] = (uint64_t)(sum6 >> 64); // 溢れたキャリーを t[7] に安全に退避
            
            // 還元（Reduction）ステップ
            uint64_t m = t[0] * n_prime;
            __uint128_t prod = (__uint128_t)m * n.limbs[0] + t[0];
            c = (uint64_t)(prod >> 64);
            
            for (int j = 1; j < 6; ++j) {
                prod = (__uint128_t)m * n.limbs[j] + t[j] + c;
                t[j - 1] = (uint64_t)prod;
                c = (uint64_t)(prod >> 64);
            }
            
            // ★ 修正3: t[7] に退避したキャリーを還元ステップの最後に合流させる
            prod = (__uint128_t)t[6] + c;
            t[5] = (uint64_t)prod;
            t[6] = (uint64_t)(prod >> 64) + t[7];
            t[7] = 0; // 次のループのためにクリア
        }
        
        // 最終的な減算
        uint384_t t_res;
        for(int i = 0; i < 6; i++) t_res.limbs[i] = t[i];
        
        if (t[6] > 0 || t_res >= n) {
            res = t_res - n;
        } else {
            res = t_res;
        }
        return res;
    }

    void add_mod(uint384_t& res, const uint384_t& a, const uint384_t& b) const {
        res = a + b;
        if (res >= n) res = res - n;
    }
    
    void sub_mod(uint384_t& res, const uint384_t& a, const uint384_t& b) const {
        uint64_t borrow = 0;
        for (int i = 0; i < 6; ++i) {
            __uint128_t diff = (__uint128_t)a.limbs[i] - b.limbs[i] - borrow;
            res.limbs[i] = (uint64_t)diff;
            borrow = (uint64_t)(diff >> 127);
        }
        if (borrow) {
            uint64_t carry = 0;
            for (int i = 0; i < 6; ++i) {
                __uint128_t sum = (__uint128_t)res.limbs[i] + n.limbs[i] + carry;
                res.limbs[i] = (uint64_t)sum;
                carry = (uint64_t)(sum >> 64);
            }
        }
    }

    void xADD(PointXZ<uint384_t>& out, const PointXZ<uint384_t>& P, const PointXZ<uint384_t>& Q, const PointXZ<uint384_t>& P_minus_Q) const {
        uint384_t t1, t2, t3, t4, t5, t6;
        sub_mod(t1, P.X, P.Z);
        add_mod(t2, Q.X, Q.Z);
        t3 = mont_mul(t1, t2);

        add_mod(t1, P.X, P.Z);
        sub_mod(t2, Q.X, Q.Z);
        t4 = mont_mul(t1, t2);

        add_mod(t5, t3, t4);
        sub_mod(t6, t3, t4);

        t1 = mont_mul(t5, t5);
        out.X = mont_mul(P_minus_Q.Z, t1);

        t2 = mont_mul(t6, t6);
        out.Z = mont_mul(P_minus_Q.X, t2);
    }

    void xDBL(PointXZ<uint384_t>& out, const PointXZ<uint384_t>& P, const uint384_t& A24_mont) const {
        uint384_t t1, t2, t3, t4, t5, t6;
        add_mod(t1, P.X, P.Z);
        sub_mod(t2, P.X, P.Z);

        t3 = mont_mul(t1, t1);
        t4 = mont_mul(t2, t2);
        sub_mod(t5, t3, t4);

        out.X = mont_mul(t3, t4);

        t6 = mont_mul(A24_mont, t5);
        add_mod(t6, t4, t6);

        out.Z = mont_mul(t5, t6);
    }
    uint384_t normal_mul(const uint384_t& a, const uint384_t& b) const {
        uint384_t temp = mont_mul(a, b);    // a * b * R^-1 mod N
        return mont_mul(temp, r2);          // (a * b * R^-1) * R^2 * R^-1 mod N = a * b mod N
    }
};

// =================================================================
// 従来の GMPベース モントゴメリ乗算コンテキスト
// =================================================================
struct MontgomeryContextGMP {
    cpp_int n;
    cpp_int r;     
    int r_bits;    
    cpp_int mask;  
    cpp_int n_prime; 
    
    // 事前確保バッファ（malloc/freeを完全に防ぐ）
    mutable mpz_t tmp_t, tmp_m, tmp_u, tmp_temp;
    mutable mpz_t t1, t2, t3, t4, t5, t6; 
    
    MontgomeryContextGMP(const cpp_int& modulus) {
        mpz_init(tmp_t); mpz_init(tmp_m); mpz_init(tmp_u); mpz_init(tmp_temp);
        mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4); mpz_init(t5); mpz_init(t6);
        
        n = modulus;
        int bits = 0;
        cpp_int temp = n;
        while (temp > 0) { bits++; temp >>= 1; }
        r_bits = ((bits + 63) / 64) * 64;
        r = cpp_int(1) << r_bits;
        mask = (cpp_int(1) << r_bits) - 1;
        
        cpp_int t = 0, newt = 1;
        cpp_int r_val = r, newr = n;
        while (newr != 0) {
            cpp_int quotient = r_val / newr;
            cpp_int temp_t = t - quotient * newt;
            t = newt; newt = temp_t;
            cpp_int temp_r = r_val - quotient * newr;
            r_val = newr; newr = temp_r;
        }
        if (t < 0) t += r;
        n_prime = (r - t) % r;
    }

    ~MontgomeryContextGMP() {
        mpz_clear(tmp_t); mpz_clear(tmp_m); mpz_clear(tmp_u); mpz_clear(tmp_temp);
        mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4); mpz_clear(t5); mpz_clear(t6);
    }

    cpp_int to_mont(const cpp_int& a) const { return (a << r_bits) % n; }
    cpp_int from_mont(const cpp_int& a) const { return mont_mul(a, 1); }

    cpp_int mont_mul(const cpp_int& a, const cpp_int& b) const {
        cpp_int res;
        mont_mul_raw(res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        return res;
    }

    void mont_mul_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_mul(tmp_t, a, b);
        mpz_mul(tmp_temp, tmp_t, n_prime.get_mpz_t());
        mpz_and(tmp_m, tmp_temp, mask.get_mpz_t());
        mpz_mul(tmp_temp, tmp_m, n.get_mpz_t());
        mpz_add(tmp_temp, tmp_t, tmp_temp);
        mpz_fdiv_q_2exp(tmp_u, tmp_temp, r_bits);
        if (mpz_cmp(tmp_u, n.get_mpz_t()) >= 0) {
            mpz_sub(tmp_u, tmp_u, n.get_mpz_t());
        }
        mpz_set(res, tmp_u);
    }

    void add_mod_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_add(res, a, b);
        if (mpz_cmp(res, n.get_mpz_t()) >= 0) mpz_sub(res, res, n.get_mpz_t());
    }

    void sub_mod_raw(mpz_t res, const mpz_t a, const mpz_t b) const {
        mpz_sub(res, a, b);
        if (mpz_sgn(res) < 0) mpz_add(res, res, n.get_mpz_t());
    }

    void xADD(PointXZ<cpp_int>& out, const PointXZ<cpp_int>& P, const PointXZ<cpp_int>& Q, const PointXZ<cpp_int>& P_minus_Q) const {
        sub_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        add_mod_raw(t2, Q.X.get_mpz_t(), Q.Z.get_mpz_t());
        mont_mul_raw(t3, t1, t2);

        add_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        sub_mod_raw(t2, Q.X.get_mpz_t(), Q.Z.get_mpz_t());
        mont_mul_raw(t4, t1, t2);

        add_mod_raw(t5, t3, t4);
        sub_mod_raw(t6, t3, t4);

        mont_mul_raw(t1, t5, t5);
        mont_mul_raw(out.X.get_mpz_t(), P_minus_Q.Z.get_mpz_t(), t1);

        mont_mul_raw(t2, t6, t6);
        mont_mul_raw(out.Z.get_mpz_t(), P_minus_Q.X.get_mpz_t(), t2);
    }

    void xDBL(PointXZ<cpp_int>& out, const PointXZ<cpp_int>& P, const cpp_int& A24_mont) const {
        add_mod_raw(t1, P.X.get_mpz_t(), P.Z.get_mpz_t());
        sub_mod_raw(t2, P.X.get_mpz_t(), P.Z.get_mpz_t());

        mont_mul_raw(t3, t1, t1);
        mont_mul_raw(t4, t2, t2);
        sub_mod_raw(t5, t3, t4);

        mont_mul_raw(out.X.get_mpz_t(), t3, t4);

        mont_mul_raw(t6, A24_mont.get_mpz_t(), t5);
        add_mod_raw(t6, t4, t6);

        mont_mul_raw(out.Z.get_mpz_t(), t5, t6);
    }

    cpp_int normal_mul(const cpp_int& a, const cpp_int& b) const {
        return mod(a * b, n);
    }
};

// --- テンプレート化: モントゴメリ・ラダーによる超高速スカラー倍 ---
// IntType と ContextType に依存し、GMP固有の処理を持たないクリーンな実装
template <typename IntType, typename ContextType>
PointXZ<IntType> mont_ladder(const PointXZ<IntType>& P, const IntType& k, const IntType& A24_mont, const ContextType& ctx, const IntType& N) {
    if (k == 0) return {ctx.to_mont(1), 0};
    
    PointXZ<IntType> R0 = P;
    PointXZ<IntType> R1;
    ctx.xDBL(R1, P, A24_mont);
    
    // ★ 追加した共通インターフェースを使用
    size_t bits = k.bit_length();
    
    for (int i = (int)bits - 2; i >= 0; --i) {
        bool bit = k.test_bit(i); // ★ 追加した共通インターフェースを使用
        if (!bit) {
            ctx.xADD(R1, R0, R1, P);    
            ctx.xDBL(R0, R0, A24_mont); 
        } else {
            ctx.xADD(R0, R0, R1, P);    
            ctx.xDBL(R1, R1, A24_mont); 
        }
    }
    return R0;
}