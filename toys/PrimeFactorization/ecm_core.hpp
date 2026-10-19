#pragma once

#include "bigint.hpp"
#include "montgomery.hpp"
#include "polynomial.hpp"
#include <emscripten/val.h>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <numeric>
//#include <iostream>

using namespace emscripten;

// =================================================================
// 状態管理構造体 (グローバルからローカルへ移行)
// =================================================================
template <typename IntType>
struct ECMState {
    IntType found_factor = IntType(0);
    long long current_prime = 0;
    IntType current_multiplier = IntType(0);
    IntType final_Z = IntType(0);
    int phase = 0;
};

// =================================================================
// 素数キャッシュとPhase 2パラメータ (C++17 inline変数で多重定義を防止)
// =================================================================
inline std::vector<int> cached_primes;
inline int max_sieved_B1 = 1;

inline int current_phase2_D = 0;
inline std::vector<int> current_phase2_coprimes;

inline void setup_phase2_params(long long B2) {
    int target_D = 210;
    if (B2 > 50000000) target_D = 510510;
    else if (B2 > 2000000) target_D = 30030;
    else if (B2 > 50000) target_D = 2310;
    
    if (current_phase2_D == target_D) return;
    
    current_phase2_D = target_D;
    current_phase2_coprimes.clear();
    
    for (int i = 1; i <= target_D / 2; i += 2) {
        // std::gcdは <numeric> のものを使用
        if (std::gcd(i, target_D) == 1) {
            current_phase2_coprimes.push_back(i);
        }
    }
}

inline void ensure_primes_up_to(int B1) {
    if (B1 <= max_sieved_B1) return;
    
    std::vector<bool> is_p(B1 + 1, true);
    is_p[0] = is_p[1] = false;
    for (int p = 2; p * p <= B1; p++) {
        if (is_p[p]) {
            for (int i = p * p; i <= B1; i += p) {
                is_p[i] = false;
            }
        }
    }
    
    cached_primes.clear();
    for (int p = 2; p <= B1; p++) {
        if (is_p[p]) cached_primes.push_back(p);
    }
    max_sieved_B1 = B1;
}

// =================================================================
// 拡張ユークリッドの互除法による逆元計算 (状態への参照を受け取る)
// =================================================================
template <typename IntType>
IntType modInverse(IntType a, IntType m, IntType& found_factor) {
    if (a == IntType(0)) {
        found_factor = m;
        return IntType(0);
    }
    
    if constexpr (std::is_same_v<IntType, cpp_int>) {
        IntType m0 = m, y = IntType(0), x = IntType(1);
        if (m == IntType(1)) return IntType(0);
        IntType q, t;
        while (a > IntType(1)) {
            if (m0 == IntType(0)) break;
            q = a / m0;
            t = m0;
            m0 = mod(a, m0);
            a = t;
            t = y;
            y = x - q * y;
            x = t;
        }
        if (a > IntType(1)) { 
            found_factor = a; 
            return IntType(1); 
        }
        if (x < IntType(0)) x += m;
        return x;
    } else {
        // ★ 符号なし固定長整数 (uint384_t) 用の安全なマイナス保護処理
        cpp_int a_gmp = to_gmp(a);
        cpp_int m_gmp = to_gmp(m);
        
        cpp_int m0 = m_gmp, y = 0, x = 1;
        if (m_gmp == 1) return IntType(0);
        cpp_int q, t;
        while (a_gmp > 1) {
            if (m0 == 0) break;
            q = a_gmp / m0;
            t = m0;
            m0 = mod(a_gmp, m0);
            a_gmp = t;
            t = y;
            y = x - q * y;
            x = t;
        }
        if (a_gmp > 1) { 
            found_factor = from_gmp(a_gmp); 
            return IntType(1); 
        }
        if (x < 0) x += m_gmp;
        return from_gmp(x);
    }
}

// =================================================================
// モントゴメリのバッチ逆元 (状態への参照渡し対応)
// =================================================================
template <typename IntType, typename ContextType>
bool batch_inverse(std::vector<IntType>& vec, const ContextType& ctx, const IntType& N, IntType& found_factor) {
    int n = vec.size();
    if (n == 0) return true;
    std::vector<IntType> prod(n);
    IntType acc = IntType(1);
    for (int i = 0; i < n; ++i) {
        prod[i] = acc;
        acc = ctx.normal_mul(acc, vec[i]); // ★ 修正
    }
    IntType acc_inv = modInverse(acc, N, found_factor);
    
    if (found_factor > IntType(0)) return false;
    
    for (int i = n - 1; i >= 0; --i) {
        IntType current = vec[i];
        vec[i] = ctx.normal_mul(acc_inv, prod[i]); // ★ 修正
        acc_inv = ctx.normal_mul(acc_inv, current); // ★ 修正
    }
    return true;
}

// =================================================================
// ECM コアエンジン (テンプレート版)
// =================================================================
template <typename IntType, typename ContextType>
val run_ecm_core(const std::string& n_str, int B1, int seed, int num_curves) {
    std::srand(seed);
    IntType N(n_str);
    val result = val::object();

    if (N % 2 == 0 || N % 3 == 0) {
        IntType factor = (N % 2 == 0) ? IntType(2) : IntType(3);
        result.set("status", "found");
        result.set("factor", factor.str());
        result.set("phase", 1);
        result.set("prime", factor.str());
        result.set("k", "1");
        result.set("Z", factor.str());
        result.set("A", "0"); result.set("B", "0");
        result.set("P_x", "0"); result.set("P_y", "0");
        result.set("curves_tried", 1);
        return result;
    }

    ensure_primes_up_to(B1);
    ContextType ctx(N); 
    
    for (int iter = 0; iter < num_curves; ++iter) {
        ECMState<IntType> state;
        state.phase = 0;
        
        // 【Suyama's Parametrization の初期化】
        IntType sigma = IntType(std::rand() % 1000000 + 6);
        IntType u = sub_mod(ctx.normal_mul(sigma, sigma), IntType(5), N);
        IntType v = ctx.normal_mul(IntType(4), sigma);
        
        IntType u2 = ctx.normal_mul(u, u);
        IntType u3 = ctx.normal_mul(u2, u);
        IntType v2 = ctx.normal_mul(v, v);
        IntType v3 = ctx.normal_mul(v2, v);
        
        PointXZ<IntType> P_XZ = {ctx.to_mont(u3), ctx.to_mont(v3)};
        
        IntType v_minus_u = sub_mod(v, u, N);
        IntType term1 = ctx.normal_mul(v_minus_u, ctx.normal_mul(v_minus_u, v_minus_u));
        IntType term2 = add_mod(ctx.normal_mul(IntType(3), u), v, N);
        IntType num = ctx.normal_mul(term1, term2);
        
        IntType den = ctx.normal_mul(IntType(16), ctx.normal_mul(u3, v));
        IntType den_inv = modInverse(den, N, state.found_factor);
        IntType A24 = ctx.normal_mul(num, den_inv);
        
        if (state.found_factor == IntType(0)) {
            IntType A24_mont = ctx.to_mont(A24);
            
            // ★ デバッグ: バッチの最初の曲線だけ、フェーズ1の開始を通知
            if (iter == 0) {
                //std::cout << "[WASM内部] 曲線 1/" << num_curves << " - フェーズ1開始 (B1=" << B1 << ", 素数個数=" << cached_primes.size() << ")" << std::endl;
            }
            
            // --- フェーズ 1 ---
            state.phase = 1; // ★ 追加: フェーズ1開始時にフラグを更新
            int prime_count = 0;
            
            for (int p : cached_primes) {
                if (p > B1) break;
                
                // ★ 超高速化 1: 巨大整数(IntType)を使わず、ネイティブの uint64_t だけで安全に上限計算を終わらせる
                uint64_t q_max = p;
                uint64_t limit = B1;
                while (q_max * p <= limit) {
                    q_max *= p;
                }
                
                P_XZ = mont_ladder(P_XZ, IntType(q_max), A24_mont, ctx, N);
                
                // ★ 超高速化 2: 毎回の重いGCD計算(GMP変換)を排除し、8192回に1回だけチェックする
                if ((++prime_count & 8191) == 0) {
                    IntType Z_normal = ctx.from_mont(P_XZ.Z);
                    IntType d = gcd(Z_normal, N);
                    if (d > IntType(1) && d < N) {
                        state.found_factor = d;
                        state.current_prime = p;
                        state.current_multiplier = IntType(q_max);
                        state.final_Z = Z_normal;
                        break;
                    } else if (d == N) {
                        state.found_factor = N;
                        break;
                    }
                }
            }

            // ★ 追加: ループを抜けたら、最後に一回だけ取りこぼしがないか確実にチェックする
            if (state.found_factor == IntType(0)) {
                IntType Z_normal = ctx.from_mont(P_XZ.Z);
                IntType d = gcd(Z_normal, N);
                if (d > IntType(1) && d < N) {
                    state.found_factor = d;
                    state.current_prime = B1; // 近似値
                    state.current_multiplier = IntType(1);
                    state.final_Z = Z_normal;
                } else if (d == N) {
                    state.found_factor = N;
                }
            }

            // --- フェーズ 2 ---
            if (state.found_factor == IntType(0)) {
                long long B2 = (long long)B1 * 500LL;
                setup_phase2_params(B2);
                
                // ★ デバッグ: フェーズ2への突入を確認
                if (iter == 0) {
                    //std::cout << "[WASM内部] フェーズ2突入 (B2=" << B2 << ", D=" << current_phase2_D << ")" << std::endl;
                }
                
                std::vector<PointXZ<IntType>> baby_steps;
                std::vector<IntType> baby_z;
                for (int j : current_phase2_coprimes) {
                    PointXZ<IntType> step = mont_ladder(P_XZ, IntType(j), A24_mont, ctx, N);
                    baby_steps.push_back(step);
                    baby_z.push_back(step.Z);
                }
                
                // ★修正: batch_inverse が見つけてもそのまま下へ流す ★ 修正: ctxを渡す
                batch_inverse(baby_z, ctx, N, state.found_factor);
                
                if (state.found_factor == IntType(0)) {
                    std::vector<IntType> baby_x;
                    for (size_t i = 0; i < baby_steps.size(); ++i) {
                        baby_x.push_back(ctx.normal_mul(baby_steps[i].X, baby_z[i])); // ★ 修正
                    }

                    Poly<IntType> F = build_poly_tree(baby_x, 0, baby_x.size() - 1, N);

                    PointXZ<IntType> DQ = mont_ladder(P_XZ, IntType(current_phase2_D), A24_mont, ctx, N);
                    
                    // ★ 修正2-A: 無限遠点 (i=0) を決して評価ツリーに混ぜない
                    long long i_min = B1 / current_phase2_D;
                    if (i_min == 0) i_min = 1; 
                    
                    long long i_max = B2 / current_phase2_D;
                    
                    if (i_min <= i_max) { // Giant-stepが存在する場合のみ実行
                        PointXZ<IntType> T_curr = mont_ladder(DQ, IntType(i_min), A24_mont, ctx, N);
                        PointXZ<IntType> T_prev;
                        if (i_min > 1) {
                            T_prev = mont_ladder(DQ, IntType(i_min - 1), A24_mont, ctx, N);
                        }
                        
                        std::vector<PointXZ<IntType>> giant_steps;
                        std::vector<IntType> giant_z;
                        
                        for (long long i = i_min; i <= i_max; ++i) {
                            giant_steps.push_back(T_curr);
                            giant_z.push_back(T_curr.Z);

                            PointXZ<IntType> T_next;
                            // ★ 修正2-B: 同じ点同士 (DQ + DQ) の場合は xADD が崩壊するため xDBL を使う
                            if (i == 1) {
                                ctx.xDBL(T_next, T_curr, A24_mont);
                            } else {
                                ctx.xADD(T_next, T_curr, DQ, T_prev);
                            }
                            T_prev = T_curr;
                            T_curr = T_next;
                        }
                        
                        // ★ 修正: ctxを渡す
                        batch_inverse(giant_z, ctx, N, state.found_factor);
                    
                        if (state.found_factor == IntType(0)) {
                            std::vector<IntType> giant_x;
                            for (size_t i = 0; i < giant_steps.size(); ++i) {
                                giant_x.push_back(ctx.normal_mul(giant_steps[i].X, giant_z[i]));
                            }

                            if (iter == 0) {
                                //std::cout << "[WASM内部] 多点評価開始: Babyサイズ=" << baby_x.size() 
                                 //         << ", Giantサイズ=" << giant_x.size() << std::endl;
                            }

                            std::vector<Poly<IntType>> eval_tree;
                            eval_tree.resize(4 * giant_x.size()); 
                            build_eval_tree(1, giant_x, 0, giant_x.size() - 1, N, eval_tree);
                            
                            // ★ 追加: ツリー構築が完走したか確認
                            //if (iter == 0) std::cout << "[WASM内部] eval_tree (部分積木) 構築完了！" << std::endl;
                            
                            std::vector<IntType> eval_results;
                            eval_results.reserve(giant_x.size()); 
                            evaluate_down(1, F, 0, giant_x.size() - 1, ctx, N, eval_tree, eval_results);

                            // ★ 追加: 割り算が完走したか確認
                            //if (iter == 0) std::cout << "[WASM内部] evaluate_down (剰余降下) 完了！" << std::endl;

                            IntType accum = IntType(1);
                            int batch_count = 0;
                            for (const IntType& res : eval_results) {
                                if (res == IntType(0)) continue;
                                accum = ctx.normal_mul(accum, res); // ★ 修正
                                
                                if (++batch_count % 128 == 0) {
                                    IntType d = gcd(accum, N);
                                    if (d > IntType(1) && d < N) {
                                        state.found_factor = d;
                                        state.phase = 2; 
                                        state.final_Z = accum;
                                        break;
                                    }
                                }
                            }
                            
                            if (state.found_factor == IntType(0) && batch_count % 128 != 0) {
                                IntType d = gcd(accum, N);
                                if (d > IntType(1) && d < N) {
                                    state.found_factor = d;
                                    state.phase = 2; 
                                    state.final_Z = accum;
                                } else if (d == N) {
                                    state.found_factor = N;
                                }
                            }
                        }
                    }
                }
            }
        } // end of protective block

        // --- 結果の返却 (スキップされずに必ず到達する) ---
        if (state.found_factor > IntType(0) && state.found_factor < N) {
            
            // ★ デバッグ: 発見の瞬間をC++側から直接出力
            //std::cout << "\n[WASM内部] ★★★ 因数ブレイクスルー発生！ ★★★" << std::endl;
            //std::cout << "  - フェーズ: " << state.phase << std::endl;
            //std::cout << "  - 因数: " << state.found_factor.str() << std::endl;
            
            IntType temp_factor = IntType(0);
            IntType A = sub_mod(ctx.normal_mul(IntType(4), A24), IntType(2), N); // ★ 修正
            IntType x = ctx.normal_mul(u3, modInverse(v3, N, temp_factor)); // ★ 修正

            result.set("status", "found");
            result.set("factor", state.found_factor.str());
            result.set("phase", state.phase);
            if (state.phase == 1) {
                result.set("prime", std::to_string(state.current_prime));
                result.set("k", state.current_multiplier.str());
            }
            result.set("Z", state.final_Z.str());
            result.set("A", A.str());
            result.set("B", "0");
            result.set("P_x", x.str());
            result.set("P_y", "0");
            result.set("curves_tried", iter + 1);
            return result; // Javascriptへ即座に返す
        }
    }
    
    result.set("status", "continue");
    result.set("curves_tried", num_curves);
    return result;
}