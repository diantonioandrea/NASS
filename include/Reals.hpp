/**
 * @file Reals.hpp
 * @author Andrea Di Antonio (github.com/diantonioandrea)
 * @brief real_t, reals_t methods.
 * @date 2024-12-11
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef NASS_REALS_HPP
#define NASS_REALS_HPP

#include <print>

#include "./Core.hpp"

#ifdef _NEON

namespace nass {
    namespace internal {

        // Load.

        /**
         * @brief Loads a reals_t from a real_t*.
         * 
         * @param Rv0 Real vector [Rv].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Ld_Rv_Rs(const real_t* Rv0) {
            #ifdef NEON16

            return vld1q_f16(Rv0);

            #else 
            #ifdef NEON32

            return vld1q_f32(Rv0);

            #else // NEON64.

            return vld1q_f64(Rv0);

            #endif
            #endif
        }

        // Store.

        /**
         * @brief Stores a reals_t into a real_t*.
         * 
         * @param Rvt0 Real vector [Rv], target [t].
         * @param Rs0 Vectorized real number [Rs].
         */
        static inline void St_RvtRs_0(real_t* Rvt0, const reals_t& Rs0) {
            #ifdef NEON16

            vst1q_f16(Rvt0, Rs0);

            #else 
            #ifdef NEON32

            vst1q_f32(Rvt0, Rs0);

            #else // NEON64.

            vst1q_f64(Rvt0, Rs0);

            #endif
            #endif
        }

        // Expand.

        /**
         * @brief Expands a real_t into a reals_t.
         * 
         * @param R0 Real number [R].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Ex_R_Rs(const real_t& R0) {
            #ifdef NEON16

            return vdupq_n_f16(R0);

            #else 
            #ifdef NEON32

            return vdupq_n_f32(R0);

            #else // NEON64.

            return vdupq_n_f64(R0);

            #endif
            #endif
        }

        // Reduce.

        /**
         * @brief Reduce a reals_t into a real_
         * 
         * @param Rs0 Real number [R].
         * @return real_t Vectorized real number [Rs].
         */
        static inline real_t Rd_Rs_R(const reals_t& Rs0) {
            #ifdef NEON16

            const float16x4_t Rs1 = vadd_f16(vget_low_f16(Rs0), vget_high_f16(Rs0));
            return vget_lane_f16(Rs1, 0) + vget_lane_f16(Rs1, 1) + vget_lane_f16(Rs1, 2) + vget_lane_f16(Rs1, 3);

            #else 
            #ifdef NEON32

            return vaddvq_f32(Rs0);

            #else // NEON64.

            return vaddvq_f64(Rs0);

            #endif
            #endif
        }

        // Operations.

        // Addition.

        /**
         * @brief Add two reals_t.
         * 
         * @param Rs0 Vectorized real number [Rs].
         * @param Rs1 Vectorized real number [Rs].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Ad_RsRs_Rs(const reals_t& Rs0, const reals_t& Rs1) {
            #ifdef NEON16

            return vaddq_f16(Rs0, Rs1);

            #else 
            #ifdef NEON32

            return vaddq_f32(Rs0, Rs1);

            #else // NEON64.

            return vaddq_f64(Rs0, Rs1);

            #endif
            #endif
        }

        // Subtraction.

        /**
         * @brief Subtracts two reals_t.
         * 
         * @param Rs0 Vectorized real number [Rs].
         * @param Rs1 Vectorized real number [Rs].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Sb_RsRs_Rs(const reals_t& Rs0, const reals_t& Rs1) {
            #ifdef NEON16

            return vsubq_f16(Rs0, Rs1);

            #else 
            #ifdef NEON32

            return vsubq_f32(Rs0, Rs1);

            #else // NEON64.

            return vsubq_f64(Rs0, Rs1);

            #endif
            #endif
        }

        // Multiplication.

        /**
         * @brief Multiplies two reals_t.
         * 
         * @param Rs0 Vectorized real number [Rs].
         * @param Rs1 Vectorized real number [Rs].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Ml_RsRs_Rs(const reals_t& Rs0, const reals_t& Rs1) {
            #ifdef NEON16

            return vmulq_f16(Rs0, Rs1);

            #else 
            #ifdef NEON32

            return vmulq_f32(Rs0, Rs1);

            #else // NEON64.

            return vmulq_f64(Rs0, Rs1);

            #endif
            #endif
        }

        // Division.

        /**
         * @brief Divides two reals_t.
         * 
         * @param Rs0 Vectorized real number [Rs].
         * @param Rs1 Vectorized real number [Rs].
         * @return reals_t Vectorized real number [Rs].
         */
        static inline reals_t Dv_RsRs_Rs(const reals_t& Rs0, const reals_t& Rs1) {
            #ifdef NEON16

            return vdivq_f16(Rs0, Rs1);

            #else 
            #ifdef NEON32

            return vdivq_f32(Rs0, Rs1);

            #else // NEON64.

            return vdivq_f64(Rs0, Rs1);

            #endif
            #endif
        }

        // Reduce (OpenMP).

        #ifdef _OPENMP
        static inline void Rd_RstRs_0(reals_t &Rst0, const reals_t &Rs0) { Rst0 = Ad_RsRs_Rs(Rst0, Rs0); }
        static inline void In_Rst_0(reals_t &Rst0) { Rst0 = Ex_R_Rs(0.0); }

        #pragma omp declare reduction(Rd_Rs: reals_t: Rd_RstRs_0(omp_out, omp_in)) initializer(In_Rst_0(omp_priv))
        #endif

    }
}

#endif

namespace nass {
    namespace internal {

        // Output.

        /**
         * @brief Prints a real_t.
         * 
         * @param R0 
         */
        static inline void Pr_R_0(const real_t& R0) {
            if(std::abs(R0) < real_tol)
                std::print(" \x1b[2m{:.3e}\033[0m ", 0.0);
            else {
                if(R0 > 0.0)
                    std::print(" ");

                std::print("{:.3e} ", R0);
            }
        }

        /**
         * @brief Prints (new line) a real_t.
         * 
         * @param R0 
         */
        static inline void Pn_R_0(const real_t& R0) {
            if(std::abs(R0) < real_tol)
                std::println(" \x1b[2m{:.3e}\033[0m ", 0.0);
            else {
                if(R0 > 0.0)
                    std::print(" ");

                std::println("{:.3e} ", R0);
            }
        }

    }
}

#endif