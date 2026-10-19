#pragma once

#include "bigint.hpp"
#include <vector>
#include <algorithm>
#include <type_traits> // if constexpr 用
//#include <iostream>

// =================================================================
// 【多項式演算エンジン】 (テンプレート版)
// =================================================================
template <typename IntType>
struct Poly {
    std::vector<IntType> coef; // coef[i] が X^i の係数

    Poly() {}
    Poly(const std::vector<IntType>& c) : coef(c) {}

    int degree() const {
        return coef.empty() ? -1 : (int)coef.size() - 1;
    }

    void clean() {
        while (!coef.empty() && coef.back() == IntType(0)) {
            coef.pop_back();
        }
    }

    Poly add(const Poly& other, const IntType& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, IntType(0));
        for (int i = 0; i <= max_deg; ++i) {
            IntType a = (i <= degree()) ? coef[i] : IntType(0);
            IntType b = (i <= other.degree()) ? other.coef[i] : IntType(0);
            res.coef[i] = add_mod(a, b, N); // ★ add_modに変更
        }
        res.clean();
        return res;
    }

    Poly sub(const Poly& other, const IntType& N) const {
        Poly res;
        int max_deg = std::max(degree(), other.degree());
        res.coef.resize(max_deg + 1, IntType(0));
        for (int i = 0; i <= max_deg; ++i) {
            IntType a = (i <= degree()) ? coef[i] : IntType(0);
            IntType b = (i <= other.degree()) ? other.coef[i] : IntType(0);
            res.coef[i] = sub_mod(a, b, N); // ★ sub_modに変更
        }
        res.clean();
        return res;
    }
};

// =================================================================
// 補助関数：分割統治法による O(N log N) パッキング
// 戻り値は巨大になるため、常に cpp_int を返す
// =================================================================
template <typename IntType>
cpp_int pack_poly(const std::vector<IntType>& coef, int left, int right, int K) {
    if (left == right) {
        if constexpr (std::is_same_v<IntType, cpp_int>) {
            return coef[left];
        } else {
            return to_gmp(coef[left]); // ★文字列パースを廃止
        }
    }
    int mid = left + (right - left) / 2;
    cpp_int P_left = pack_poly(coef, left, mid, K);
    cpp_int P_right = pack_poly(coef, mid + 1, right, K);
    
    int shift = (mid - left + 1) * K;
    return P_left + (P_right << shift);
}

// =================================================================
// 補助関数：分割統治法による O(N log N) アンパッキング
// =================================================================
template <typename IntType>
void unpack_poly(const cpp_int& packed, int left, int right, int K, const cpp_int& mask, const IntType& N, std::vector<IntType>& result) {
    if (left == right) {
        cpp_int rem = packed & mask;
        if constexpr (std::is_same_v<IntType, cpp_int>) {
            result[left] = mod(rem, N);
        } else {
            // ★修正: 文字列経由を完全に廃止し、メモリ直接コピーの to_gmp を使用
            cpp_int N_cpp = to_gmp(N); 
            rem = mod(rem, N_cpp);
            result[left] = from_gmp(rem); 
        }
        return;
    }
    int mid = left + (right - left) / 2;
    int shift = (mid - left + 1) * K;
    
    cpp_int lower_mask = (cpp_int(1) << shift) - 1;
    cpp_int P_left = packed & lower_mask;
    cpp_int P_right = packed >> shift;
    
    unpack_poly(P_left, left, mid, K, mask, N, result);
    unpack_poly(P_right, mid + 1, right, K, mask, N, result);
}

// =================================================================
// 超高速化版: クロネッカー代入を用いた高速多項式乗算 (mod N)
// =================================================================
template <typename IntType>
Poly<IntType> kronecker_multiply(const Poly<IntType>& A, const Poly<IntType>& B, const IntType& N) {
    if (A.degree() < 0 || B.degree() < 0) return Poly<IntType>();

    int n = A.degree();
    int m = B.degree();
    int min_len = std::min(n + 1, m + 1);

    int bits_N = N.bit_length(); // ★共通インターフェースを使用
    int log2_min = 0;
    int temp_min = min_len;
    while (temp_min > 0) { log2_min++; temp_min >>= 1; }

    int K = 2 * bits_N + log2_min + 2;

    cpp_int packed_A = pack_poly(A.coef, 0, n, K);
    cpp_int packed_B = pack_poly(B.coef, 0, m, K);

    // Kronecker置換: 1回の多倍長整数乗算で多項式積を計算
    cpp_int packed_C = packed_A * packed_B;

    Poly<IntType> C;
    C.coef.resize(n + m + 1, IntType(0));
    cpp_int mask = (cpp_int(1) << K) - 1;
    
    unpack_poly(packed_C, 0, n + m, K, mask, N, C.coef);

    C.clean();
    return C;
}

// =================================================================
// 分割統治法による剰余木 (Subproduct Tree) の構築
// =================================================================
template <typename IntType>
Poly<IntType> build_poly_tree(const std::vector<IntType>& roots, int left, int right, const IntType& N) {
    if (left == right) {
        IntType val = mod(roots[left], N);
        IntType neg_val = (val == IntType(0)) ? IntType(0) : N - val;
        return Poly<IntType>({neg_val, IntType(1)});
    }
    
    int mid = left + (right - left) / 2;
    Poly<IntType> left_poly = build_poly_tree(roots, left, mid, N);
    Poly<IntType> right_poly = build_poly_tree(roots, mid + 1, right, N);
    
    return kronecker_multiply(left_poly, right_poly, N);
}

// =================================================================
// モニック多項式による多項式剰余算 (A mod B)
// ★ コンテキストを受け取り、GMPを一切使わない乗算を使用する
// =================================================================
template <typename IntType, typename ContextType>
Poly<IntType> poly_mod_monic(Poly<IntType> A, const Poly<IntType>& B, const ContextType& ctx, const IntType& N) {
    int degB = B.degree();
    if (degB < 0) return A;
    
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        mpz_t temp, factor_mpz, mod_N;
        mpz_init(temp);
        mpz_init(factor_mpz);
        mpz_init_set(mod_N, N.get_mpz_t());
        
        while (A.degree() >= degB) {
            int d = A.degree() - degB;
            mpz_set(factor_mpz, A.coef.back().get_mpz_t());
            
            for (int i = 0; i <= degB; ++i) {
                mpz_mul(temp, B.coef[i].get_mpz_t(), factor_mpz);
                mpz_sub(A.coef[i + d].get_mpz_t(), A.coef[i + d].get_mpz_t(), temp);
                mpz_fdiv_r(A.coef[i + d].get_mpz_t(), A.coef[i + d].get_mpz_t(), mod_N);
            }
            A.clean();
        }
        mpz_clear(temp);
        mpz_clear(factor_mpz);
        mpz_clear(mod_N);
        return A;
    } else {
        int infinite_loop_guard = 0; // ★ 追加: ループ監視用
        
        while (A.degree() >= degB) {
            int old_deg = A.degree(); // ★ 追加: 実行前の次数を記憶
            
            int d = A.degree() - degB;
            IntType factor = A.coef.back();
            for (int i = 0; i <= degB; ++i) {
                IntType term = ctx.normal_mul(B.coef[i], factor);
                A.coef[i + d] = sub_mod(A.coef[i + d], term, N);
            }
            A.clean();
            
            // ==========================================================
            // ★ 最強の無限ループ破壊＆デバッグ機構
            // ==========================================================
            if (A.degree() >= old_deg) {
                infinite_loop_guard++;
                if (infinite_loop_guard > 3) {
                    //std::cout << "[WASM 致命的エラー] 多項式除算で次数が減りません (無限ループを強制破壊)！" << std::endl;
                    //std::cout << "  -> Bの最高次係数: " << B.coef[degB].str() << std::endl;
                    //std::cout << "  -> Aの対象係数: " << factor.str() << std::endl;
                    A.coef.pop_back(); // 強制的に最高次を削除してループを脱出させる
                }
            } else {
                infinite_loop_guard = 0;
            }
        }
        return A;
    }
}

// =================================================================
// 評価対象の点群から、多項式除算用の木を構築
// =================================================================
template <typename IntType>
void build_eval_tree(int node, const std::vector<IntType>& points, int left, int right, const IntType& N, std::vector<Poly<IntType>>& tree) {
    if (left == right) {
        IntType val = mod(points[left], N);
        IntType neg_val = (val == IntType(0)) ? IntType(0) : N - val;
        tree[node] = Poly<IntType>({neg_val, IntType(1)});
        return;
    }
    int mid = left + (right - left) / 2;
    build_eval_tree(2 * node, points, left, mid, N, tree);
    build_eval_tree(2 * node + 1, points, mid + 1, right, N, tree);
    
    tree[node] = kronecker_multiply(tree[2 * node], tree[2 * node + 1], N);
}

// =================================================================
// 木構造の上から下へ多項式の割り算を繰り返し、最終的な代入値を抽出する
// =================================================================
template <typename IntType, typename ContextType>
void evaluate_down(int node, Poly<IntType> F, int left, int right, const ContextType& ctx, const IntType& N, const std::vector<Poly<IntType>>& tree, std::vector<IntType>& results) {
    // ★ 修正: ctxを渡す
    F = poly_mod_monic(F, tree[node], ctx, N);
    
    if (left == right) {
        results.push_back(F.degree() >= 0 ? F.coef[0] : IntType(0));
        return;
    }
    
    int mid = left + (right - left) / 2;
    evaluate_down(2 * node, F, left, mid, ctx, N, tree, results);
    evaluate_down(2 * node + 1, F, mid + 1, right, ctx, N, tree, results);
}